import errno
import os
import shutil
from subprocess import Popen, PIPE

import pandas as pd


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

    return mapping.get(extension, None)


def _get_all_files(dirpath):
    files = []
    for root, dirnames, filenames in os.walk(dirpath):
        for filename in filenames:
            files.append(os.path.join(root, filename))
    return files


def _get_new_rgid(rgid, old_sample_id, new_sample_id):
    assert len(rgid.split('_')) == 4
    sample, library, lane, flowcell = rgid.split('_')
    assert sample.replace('_', '').replace('-', '') == old_sample_id or 'DLPNegative' in sample or 'DLPGm' in sample
    rgid = '_'.join((new_sample_id, library, lane, flowcell))
    return rgid


def _update_rg_line(rg, old_sample_id, new_sample_id):
    rg = rg.strip().split()

    assert rg[0] == '@RG'
    assert rg[1].startswith('ID:')

    rgid = rg[1][len('ID:'):]
    rg[1] = 'ID:' + _get_new_rgid(rgid, old_sample_id, new_sample_id)

    sample_idx = [i for i, v in enumerate(rg) if v.startswith('SM:')]
    assert len(sample_idx) == 1
    sample_idx = sample_idx[0]

    assert rg[sample_idx][len('SM:'):] == old_sample_id

    rg[sample_idx] = 'SM:{}'.format(new_sample_id)

    return '\t'.join(rg) + '\n'


def _update_barcode_comment(barcode_line):
    assert barcode_line.startswith('@CO\tCB:')

    barcode = barcode_line[len('@CO\tCB:'):]
    assert len(barcode.split('-')) == 4
    barcode = '-'.join(barcode.split('-')[1:])

    return '@CO\tCB:{}\n'.format(barcode)


def _update_read(read, old_sample_id, new_sample_id):
    read = read.strip().split('\t')

    barcode_idx = [i for i, v in enumerate(read) if v.startswith('CB:Z:')]
    assert len(barcode_idx) == 1
    barcode_idx = barcode_idx[0]
    read[barcode_idx] = 'CB:Z:' + '-'.join(read[barcode_idx][len('CB:Z:'):].split('-')[1:])

    rg_idx = [i for i, v in enumerate(read) if v.startswith('RG:Z:')]
    assert len(rg_idx) == 1
    rg_idx = rg_idx[0]

    read[rg_idx] = 'CB:Z:' + _get_new_rgid(read[rg_idx], old_sample_id, new_sample_id)

    return '\t'.join(read) + '\n'


def update_bam(bamfile, tempdir, old_sample_id, new_sample_id, output_bam):
    if os.path.getsize(bamfile) == 7:
        print('empty bam: {}'.format(bamfile))
        shutil.copyfile(bamfile, output_bam)
        return

    sampath = os.path.join(tempdir, 'converted.sam')
    cmd = ['samtools', 'view', '-h', '-o', sampath, bamfile]
    _run_cmd(cmd)

    updated_sam = os.path.join(tempdir, 'updated.sam')
    with open(sampath, 'rt') as reader, open(updated_sam, 'wt') as writer:
        for line in reader:

            if not line.startswith('@'):
                line = _update_read(line, old_sample_id, new_sample_id)
            elif line.startswith('@RG'):
                line = _update_rg_line(line, old_sample_id, new_sample_id)
            elif line.startswith('@CO'):
                line = _update_barcode_comment(line)

            writer.write(line)

    cmd = ['samtools', 'view', '-bSh', '-o', output_bam, updated_sam]
    _run_cmd(cmd)

    cmd = ['samtools', 'index', output_bam]
    _run_cmd(cmd)


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
        return '-'.join(cell_id.split('-')[1:])


def update_csv(filepath, old_sample_id, new_sample_id, output_csv):
    df = pd.read_csv(filepath)

    if 'sample_id' in df.columns.values:
        df['sample_id'] = df['sample_id'].apply(lambda x: _update_sample(x, old_sample_id, new_sample_id))

    if 'cell_id' in df.columns.values:
        df['cell_id'] = df['cell_id'].apply(lambda x: _update_cell_id(x))

    df.to_csv(output_csv, index=False)


def mondrian(old_sample_id, new_sample_id, input_dir, output_dir, tempdir):
    _makedirs(tempdir)
    _makedirs(output_dir)

    allfiles = _get_all_files(input_dir)

    for filepath in allfiles:
        output_path = os.path.join(output_dir, os.path.basename(filepath))
        if _get_file_type(filepath) == 'bam':
            update_bam(filepath, tempdir, old_sample_id, new_sample_id, output_path)
        elif _get_file_type(filepath) == 'csv':
            update_csv(filepath, old_sample_id, new_sample_id, output_path)
        elif _get_file_type(filepath) in ['csv_yaml', 'pdf', 'tar']:
            shutil.copyfile(filepath, output_path)
        else:
            print('skipping: {}'.format(filepath))


mondrian('RPE-Noco', 'A12345', 'results', 'results_updated', 'tempdir')
