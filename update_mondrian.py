import errno
import os
import shutil
import sys
from subprocess import Popen, PIPE

import pandas as pd
import pysam
import tqdm
import yaml


def _load_flagstat(flagstat_file):
    data = {}
    with open(flagstat_file, 'rt') as reader:
        for line in reader:
            line_split = line.strip().split()

            if 'total (QC-passed reads + QC-failed reads)' in line:
                data['total'] = (line_split[0], line_split[2])
            elif 'secondary' in line:
                data['secondary'] = (line_split[0], line_split[2])
            elif 'supplementary' in line:
                data['supplementary'] = (line_split[0], line_split[2])
            elif 'duplicates' in line:
                data['duplicates'] = (line_split[0], line_split[2])
            elif 'mapped' in line:
                data['mapped'] = (line_split[0], line_split[2])
            elif 'paired in sequencing' in line:
                data['paired'] = (line_split[0], line_split[2])
            elif 'read1' in line:
                data['read1'] = (line_split[0], line_split[2])
            elif 'read2' in line:
                data['read2'] = (line_split[0], line_split[2])
            elif 'properly paired' in line:
                data['properly_paired'] = (line_split[0], line_split[2])
            elif 'singletons' in line:
                data['singletons'] = (line_split[0], line_split[2])
    return data


def compare_bams(oldbam, newbam, tempdir):
    old_fstat = os.path.join(tempdir, 'old_flagstat.txt')
    new_fstat = os.path.join(tempdir, 'new_flagstat.txt')

    cmd = ['samtools', 'flagstat', oldbam]
    _run_cmd(cmd, output=old_fstat)

    cmd = ['samtools', 'flagstat', newbam]
    _run_cmd(cmd, output=new_fstat)

    old_fstat = _load_flagstat(old_fstat)
    new_fstat = _load_flagstat(new_fstat)

    assert old_fstat == new_fstat


def _makedirs(dirname):
    dirname = os.path.abspath(dirname)
    try:
        os.makedirs(dirname)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    assert os.path.isdir(dirname)


def _run_cmd(cmd, output=None):
    stdout = None
    if output:
        stdout = open(output, "w")

    p = Popen(cmd, stdout=stdout, stderr=PIPE)

    cmdout, cmderr = p.communicate()

    retc = p.returncode

    if retc:
        raise Exception(
            "command failed. stderr:{}, stdout:{}".format(cmdout, cmderr)
        )

    if output:
        stdout.close()


def _get_file_type(filepath):
    extension = filepath[filepath.index('.'):]

    mapping = {
        '.csv.gz': 'csv',
        '.csv.gz.yaml': 'csv_yaml',
        '.bam': 'bam',
        '.pdf': 'pdf',
        '.tar': 'tar',
        '.tar.gz': 'tar'
    }

    if os.path.basename(filepath) == 'metadata.yaml':
        return 'meta_yaml'

    return mapping.get(extension, None)


def _get_all_files(dirpath):
    files = []
    for root, dirnames, filenames in os.walk(dirpath):
        for filename in filenames:
            files.append(os.path.join(root, filename))
    return files


def _update_bam_header(header, old_sample_id, new_sample_id):
    # allowed_samples = {'DLPNegative', 'DLPGm', old_sample_id}

    comments = header['CO']
    updated_comments = []
    for comment in comments:
        if comment.startswith('CB:'):
            comment = comment[len('CB:'):]
            comment = comment.rsplit('-', maxsplit=3)
            comment = '-'.join(comment[1:])
            comment = 'CB:' + comment
        updated_comments.append(comment)
    header['CO'] = updated_comments

    rg_rename = {}

    for rg in header['RG']:
        sample, library, lane, flowcell = rg['ID'].rsplit('_', maxsplit=3)

        # assert sample in allowed_samples
        assert rg['ID'] not in rg_rename

        rg_id = f'{new_sample_id}_{library}_{lane}_{flowcell}'

        rg_rename[rg['ID']] = rg_id

        rg['ID'] = rg_id
        rg['PU'] = f'{lane}_{flowcell}'
        rg['SM'] = new_sample_id

    return rg_rename


def _update_bam_reads(input_bam, output_bam, old_sample_id, rg_mapping):
    # allowed_samples = {'DLPNegative', 'DLPGm', old_sample_id}

    for a in tqdm.tqdm(input_bam, total=input_bam.mapped + input_bam.unmapped):
        tags = []
        for t in a.tags:
            if t[0] == 'RG':
                tags.append(('RG', rg_mapping[t[1]]))
            elif t[0] == 'CB':
                cell_id = t[1]
                cell_info = cell_id.rsplit('-', maxsplit=3)
                # assert cell_info[0] in allowed_samples
                new_cell_id = '-'.join(cell_info[1:])
                tags.append(('CB', new_cell_id))
            else:
                tags.append(t)
        a.tags = tags
        output_bam.write(a)


def update_bam(bamfile, output_bam, old_sample_id, new_sample_id, tempdir):
    if os.path.getsize(bamfile) == 7:
        print('empty bam: {}'.format(bamfile))
        shutil.copyfile(bamfile, output_bam)
        shutil.copyfile(bamfile + '.bai', output_bam + '.bai')
        return

    bam = pysam.AlignmentFile(bamfile, 'rb')
    header = bam.header.to_dict()

    rg_mapping = _update_bam_header(header, old_sample_id, new_sample_id)

    assert len([v['ID'] for v in header['RG']]) == len(set([v['ID'] for v in header['RG']]))

    with pysam.AlignmentFile(output_bam, 'wb', header=header) as output:
        _update_bam_reads(bam, output, old_sample_id, rg_mapping)

    bam.close()

    cmd = ['samtools', 'index', output_bam]
    _run_cmd(cmd)

    compare_bams(bamfile, output_bam, tempdir)


def _update_sample(curr_sample_id, old_sample_id, new_sample_id):
    if 'DLPNegative' in curr_sample_id or 'DLPGm' in curr_sample_id:
        return curr_sample_id
    else:
        assert curr_sample_id == old_sample_id, (curr_sample_id, old_sample_id)
        return new_sample_id


def _update_cell_id(cell_id):
    if cell_id == 'reference':
        return cell_id
    else:
        return '-'.join(cell_id.rsplit('-', maxsplit=3)[1:])


def update_csv(filepath, old_sample_id, new_sample_id, output_csv):
    df = pd.read_csv(filepath)

    if 'sample_id' in df.columns.values:
        samples = list(df[~ df['is_control']]['sample_id'].unique())
        assert len(samples) == 1
        assert samples[0] == old_sample_id
        df['sample_id'] = new_sample_id

    if 'cell_id' in df.columns.values:
        df['cell_id'] = df['cell_id'].apply(lambda x: _update_cell_id(x))

    df.to_csv(output_csv, index=False)


def update_metadata_yaml(meta_yaml_input, meta_yaml_output, old_sample, new_sample):
    with open(meta_yaml_input, 'rt') as reader:
        data = yaml.safe_load(reader)

    cells = data['meta']['cell_ids']
    cells = ['-'.join(v.rsplit('-', maxsplit=3)[1:]) for v in cells]
    data['meta']['cell_ids'] = cells

    assert old_sample in data['meta']['sample_ids']
    data['meta']['sample_ids'] = [new_sample]

    with open(meta_yaml_output, 'wt') as writer:
        yaml.dump(data, writer, default_flow_style=False)


def mondrian(old_sample_id, new_sample_id, input_dir, output_dir, tempdir):
    _makedirs(tempdir)
    _makedirs(output_dir)

    allfiles = _get_all_files(input_dir)

    for filepath in allfiles:

        output_path = os.path.join(output_dir, os.path.basename(filepath))
        if _get_file_type(filepath) == 'bam':
            update_bam(filepath, output_path, old_sample_id, new_sample_id, tempdir)
        elif _get_file_type(filepath) == 'csv':
            update_csv(filepath, old_sample_id, new_sample_id, output_path)
        elif _get_file_type(filepath) in ['csv_yaml', 'pdf', 'tar']:
            shutil.copyfile(filepath, output_path)
        elif _get_file_type(filepath) == 'meta_yaml':
            update_metadata_yaml(filepath, output_path, old_sample_id, new_sample_id)
        else:
            print('skipping: {}'.format(filepath))


old_sample_id = sys.argv[1]
new_sample_id = sys.argv[2]
inputdir = sys.argv[3]
outputdir = sys.argv[4]
tempdir = sys.argv[5]

mondrian(
    old_sample_id, new_sample_id, inputdir, outputdir, tempdir
)

# mondrian(
#     'RPE-Noco', 'A12345',
#     '/juno/work/shah/isabl_data_lake/analyses/40/15/24015/results/',
#     'results_updated', 'tempdir'
# )
