import errno
import gzip
import os
import shutil

import pandas as pd


def _get_separator(filepath):
    opener = gzip.open if filepath.endswith('.gz') else open

    with opener(filepath, 'rt') as reader:
        header = reader.readline()
        # print(header)

        # if ',' in header and '\t' in header:
        #     raise Exception()

        if '\t' in header:
            return '\t'
        elif ',' in header:
            return ','
        else:
            raise Exception('unknown separator in {}'.format(filepath))


def _makedirs(dirname):
    dirname = os.path.abspath(dirname)
    try:
        os.makedirs(dirname)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
    assert os.path.isdir(dirname)


def _run_cmd(cmd):
    cmd = ' '.join(cmd)
    os.system(cmd)


def _get_file_type(filepath):
    extension = filepath[filepath.index('.'):]

    mapping = {
        '.csv.gz': 'csv',
        '.csv.gz.yaml': 'csv_yaml',
        '.bam': 'bam',
        '.pdf': 'pdf',
        '.maf': 'maf',
        '.tar': 'tar',
        '.tar.gz': 'tar',
        '.vcf': 'vcf',
        '.vcf.gz': 'vcf',
        '.seg': 'csv'
    }

    return mapping.get(extension, None)


def _get_all_files(dirpath):
    files = os.listdir(dirpath)
    files = [os.path.join(dirpath, v) for v in files]
    return files


def _update_rg_line(rg, old_sample_id, new_sample_id):
    rg = rg.strip().split()

    assert rg[0] == '@RG'
    assert rg[1].startswith('ID:')

    sample_idx = [i for i, v in enumerate(rg) if v.startswith('SM:')]
    assert len(sample_idx) == 1
    sample_idx = sample_idx[0]

    assert rg[sample_idx][len('SM:'):] == old_sample_id

    rg[sample_idx] = 'SM:{}'.format(new_sample_id)

    return '\t'.join(rg) + '\n'


def update_bam(bamfile, tempdir, old_sample_id, new_sample_id, output_bam):
    headerfile = os.path.join(tempdir, 'header.sam')
    cmd = ['samtools', 'view', '-H', '-o', headerfile, bamfile]
    _run_cmd(cmd)

    updated_header = os.path.join(tempdir, 'updated_header.sam')
    with open(headerfile, 'rt') as reader, open(updated_header, 'wt') as writer:
        for line in reader:
            if line.startswith('@RG'):
                line = _update_rg_line(line, old_sample_id, new_sample_id)

            writer.write(line)

    cmd = ['samtools', 'reheader', updated_header, bamfile, '>', output_bam]
    _run_cmd(cmd)

    cmd = ['samtools', 'index', output_bam]
    _run_cmd(cmd)


def update_csv(filepath, old_sample_id, new_sample_id, output_csv, sample_col='sample_id'):
    sep = _get_separator(filepath)

    df = pd.read_csv(filepath, sep=sep)
    assert df[sample_col].iloc[0] == old_sample_id
    df[sample_col] = new_sample_id
    df.to_csv(output_csv, sep=sep, index=False)


def update_germline_vcf(filepath, old_sample_id, new_sample_id, output_path, tempdir):
    opener = gzip.open if filepath.endswith('.gz') else open

    temp_output_path = os.path.join(tempdir, 'temp_germline_output.vcf')

    with opener(filepath, 'rt') as reader, open(temp_output_path, 'wt') as writer:
        for line in reader:
            if line.startswith('#') and not line.startswith('##'):
                line = line.strip().split('\t')
                assert line[-1] == old_sample_id, (line, old_sample_id)
                line[-1] = new_sample_id
                line = '\t'.join(line) + '\n'
            writer.write(line)

    if filepath.endswith('.gz'):
        cmd = ['bgzip', temp_output_path, '-c', '>', output_path]
        _run_cmd(cmd)
        cmd = ['tabix', '-f', '-p', 'vcf', output_path]
        _run_cmd(cmd)
        cmd = ['bcftools', 'index', output_path]
        _run_cmd(cmd)
    else:
        shutil.copyfile(temp_output_path, output_path)


def update_paired_vcf(filepath, old_normal_id, new_normal_id, old_tumour_id, new_tumour_id, output_path, tempdir):
    opener = gzip.open if filepath.endswith('.gz') else open

    temp_output_path = os.path.join(tempdir, 'temp_germline_output.vcf')

    with opener(filepath, 'rt') as reader, open(temp_output_path, 'wt') as writer:
        for line in reader:
            if line.startswith('#') and not line.startswith('##'):
                line = line.replace(old_tumour_id, new_tumour_id)
                line = line.replace(old_normal_id, new_normal_id)
            writer.write(line)

    if filepath.endswith('.gz'):
        cmd = ['bgzip', temp_output_path, '-c', '>', output_path]
        _run_cmd(cmd)
        cmd = ['tabix', '-f', '-p', 'vcf', output_path]
        _run_cmd(cmd)
        cmd = ['bcftools', 'index', output_path]
        _run_cmd(cmd)
    else:
        shutil.copyfile(temp_output_path, output_path)


def update_maf(filepath, old_normal_id, new_normal_id, old_tumour_id, new_tumour_id, output_path, tempdir):
    with open(filepath) as infile_read:
        maf_header = infile_read.readline()
    assert maf_header.startswith('#version 2.4')

    df = pd.read_csv(filepath, dtype='str', skiprows=1, sep='\t')

    assert len(df['Tumor_Sample_Barcode'].unique()) == 1
    assert df['Tumor_Sample_Barcode'].unique()[0] == old_tumour_id
    df['Tumor_Sample_Barcode'] = new_tumour_id

    assert len(df['Matched_Norm_Sample_Barcode'].unique()) == 1
    assert df['Matched_Norm_Sample_Barcode'].unique()[0] == old_normal_id
    df['Matched_Norm_Sample_Barcode'] = new_normal_id

    with open(output_path, 'wt') as writer:
        writer.write('#version 2.4\n')

        df.to_csv(writer, sep='\t', index=False)


def alignment(old_sample_id, new_sample_id, input_dir, output_dir, tempdir):
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


def germline(old_sample_id, new_sample_id, input_dir, output_dir, tempdir):
    _makedirs(tempdir)
    _makedirs(output_dir)

    allfiles = _get_all_files(input_dir)

    for filepath in allfiles:
        output_path = os.path.join(output_dir, os.path.basename(filepath))

        if _get_file_type(filepath) == 'vcf':
            update_germline_vcf(filepath, old_sample_id, new_sample_id, output_path, tempdir)
        elif _get_file_type(filepath) == 'csv':
            update_csv(filepath, old_sample_id, new_sample_id, output_path, sample_col='sample')
        elif _get_file_type(filepath) in ['csv_yaml', 'pdf', 'tar', 'maf']:
            shutil.copyfile(filepath, output_path)
        else:
            print('skipping: {}'.format(filepath))


def breakpoint(old_tumour_id, new_tumour_id, old_normal_id, new_normal_id, input_dir, output_dir, tempdir):
    _makedirs(tempdir)
    _makedirs(output_dir)

    allfiles = _get_all_files(input_dir)

    for filepath in allfiles:
        output_path = os.path.join(output_dir, os.path.basename(filepath))

        if _get_file_type(filepath) == 'vcf':
            update_paired_vcf(
                filepath, old_normal_id, new_normal_id, old_tumour_id, new_tumour_id, output_path,
                tempdir
            )
        elif _get_file_type(filepath) in ['csv_yaml', 'pdf', 'tar', 'maf', 'csv']:
            shutil.copyfile(filepath, output_path)
        else:
            print('skipping: {}'.format(filepath))


def somatic(old_tumour_id, new_tumour_id, old_normal_id, new_normal_id, input_dir, output_dir, tempdir):
    _makedirs(tempdir)
    _makedirs(output_dir)

    allfiles = _get_all_files(input_dir)

    for filepath in allfiles:
        output_path = os.path.join(output_dir, os.path.basename(filepath))

        if _get_file_type(filepath) == 'vcf':
            update_paired_vcf(
                filepath, old_normal_id, new_normal_id, old_tumour_id, new_tumour_id, output_path,
                tempdir
            )
        if _get_file_type(filepath) == 'maf':
            update_maf(
                filepath, old_normal_id, new_normal_id, old_tumour_id, new_tumour_id, output_path,
                tempdir
            )
        elif _get_file_type(filepath) in ['csv_yaml', 'pdf', 'tar', 'maf', 'csv']:
            shutil.copyfile(filepath, output_path)
        else:
            print('skipping: {}'.format(filepath))


def titan(old_tumour_id, new_tumour_id, old_normal_id, new_normal_id, input_dir, output_dir, tempdir):
    _makedirs(tempdir)
    _makedirs(output_dir)

    allfiles = _get_all_files(input_dir)

    for filepath in allfiles:
        output_path = os.path.join(output_dir, os.path.basename(filepath))

        if _get_file_type(filepath) == 'vcf':
            update_paired_vcf(
                filepath, old_normal_id, new_normal_id, old_tumour_id, new_tumour_id, output_path,
                tempdir
            )
        elif _get_file_type(filepath) == 'csv':
            if filepath.endswith('.seg'):
                update_csv(filepath, old_tumour_id, new_tumour_id, output_path, sample_col='sample')
            elif 'titan_parsed' in filepath or 'titan_segs' in filepath:
                update_csv(filepath, old_tumour_id, new_tumour_id, output_path, sample_col='Sample')
            else:
                shutil.copyfile(filepath, output_path)
        elif _get_file_type(filepath) in ['csv_yaml', 'pdf', 'tar']:
            shutil.copyfile(filepath, output_path)
        else:
            print('skip:{}'.format(filepath))


alignment(
    'AK-RT-020-N', 'A12345',
    '/juno/work/shah/isabl_data_lake/analyses/32/74/23274/results/AK-RT-020-N_NORMAL/', 'results_updated',
    'tempdir'
)
# germline(
#     'AK-RT-020-N', 'A12345',
#     '/juno/work/shah/isabl_data_lake/analyses/38/35/23835/results/germline/AK-RT-020-N/', 'results_updated',
#     'tempdir'
# )
# breakpoint(
#     'AK-RT-020-T', 'T12345', 'AK-RT-020-N', 'N12345',
#     '/juno/work/shah/isabl_data_lake/analyses/45/73/24573/results/breakpoints/AK-RT-020-T', 'results_updated',
#     'tempdir'
# )
# somatic(
#     'AK-RT-020-T', 'T12345', 'AK-RT-020-N', 'N12345',
#     '/juno/work/shah/isabl_data_lake/analyses/37/84/23784/results/somatic/AK-RT-020-T', 'results_updated',
#     'tempdir'
# )
# titan(
#     'AK-RT-020-T', 'T12345', 'AK-RT-020-N', 'N12345',
#     '/juno/work/shah/isabl_data_lake/analyses/37/96/23796/results/copynumber/AK-RT-020-T/titan', 'results_updated',
#     'tempdir'
# )
