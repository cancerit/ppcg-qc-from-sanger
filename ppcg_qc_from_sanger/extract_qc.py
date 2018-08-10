import os
import shutil
import sys
import errno
import tarfile
import logging
import json
import re
import gzip
from tempfile import TemporaryDirectory
from . import exec_subp_and_wait


# setup logs
logger = logging.getLogger('ppcg-qc-from-sanger-extract_qc')
# create console handler and set level to debug
cha = logging.StreamHandler()
cha.setLevel(logging.DEBUG)
cha.setFormatter(logging.Formatter('%(asctime)s %(levelname)8s - %(message)s',
                                   '%Y-%m-%d %H:%M:%S'))
# add ch to logger
logger.addHandler(cha)
logger.setLevel(logging.DEBUG)

BAS_HEADER = [
    'bam_filename',
    'sample',
    'platform',
    'platform_unit',
    'library',
    'readgroup',
    'read_length_r1',
    'read_length_r2',
    '#_mapped_bases',
    '#_mapped_bases_r1',
    '#_mapped_bases_r2',
    '#_divergent_bases',
    '#_divergent_bases_r1',
    '#_divergent_bases_r2',
    '#_total_reads',
    '#_total_reads_r1',
    '#_total_reads_r2',
    '#_mapped_reads',
    '#_mapped_reads_r1',
    '#_mapped_reads_r2',
    '#_mapped_reads_properly_paired',
    '#_gc_bases_r1',
    '#_gc_bases_r2',
    'mean_insert_size',
    'insert_size_sd',
    'median_insert_size',
    '#_duplicate_reads',
    '#_mapped_pairs',
    '#_inter_chr_pairs']

OUTPUT_HEADER = [
    'Sample name',
    'Sample Type',
    'Depth',
    'Fraction of mapped reads',
    'Insert size',
    'Insert size sd',
    'r1 GC content',
    'r2 GC content',
    'Fraction of duplicated reads',
    'Fraction of mis-matched pairs',
    'Contamination',
    'Gender',
    'Fraction of matched gender',
    'Fraction of matched genotype',
    'Normal contamination']

VARIANT_COUNT_HEADER = [
    'Number of SNVs', 'Number of INDELs', 'Number of SVs', 'Number of CNVs']

# global variable to determine if need to count variants
count_variants = False


def extract_from_sanger(args):
    '''
    the main function for handling the whole QC metrics extraction process
    '''
    if not args.debug:
        logger.setLevel(logging.INFO)

    tumour_bas = get_abs_path(args.tumour_bas)
    normal_bas = get_abs_path(args.normal_bas)
    genome_size = args.genome_size
    if not isinstance(genome_size, int):
        raise RuntimeError('genome_size is not int')
    variant_call_tar = get_abs_path(args.variant_call_tar)
    output_tar = get_abs_path(args.output_tar)

    # chcck files existence
    check_file_exists(tumour_bas)
    check_file_exists(normal_bas)
    check_file_exists(variant_call_tar)

    # test if output is writable
    if not os.path.exists(output_tar):
        try:
            with open(output_tar, 'w') as out:
                out.write('place holder\n')
        except OSError as exc:
            raise RuntimeError('output is not writable: %s.' % str(exc))
        finally:
            os.remove(output_tar)
    else:
        raise RuntimeError('existing output file: %s.' % output_tar)

    # if bas files hava valid columns
    validate_bas(tumour_bas)
    validate_bas(normal_bas)
    validate_tar_name(variant_call_tar)
    validate_tar_name(output_tar)

    t_bas_content = get_bas_content(tumour_bas)
    t_sample_name = get_sample_name_from_bas(t_bas_content)
    n_bas_content = get_bas_content(normal_bas)
    n_sample_name = get_sample_name_from_bas(n_bas_content)

    # determine if needs to count variants
    global count_variants
    count_variants = args.count_variants

    # create a temp dir to store extracted files
    with TemporaryDirectory() as temp_dir:
        # extracted required files and move genotyping files to output folder
        (extracted_gender_file,
         extracted_purity_file,
         extracted_contamination_files,
         extracted_genotype_files,
         extracted_snv_file,
         extracted_indel_file,
         extracted_sv_file,
         extracted_cnv_file) = \
            extract_and_place_required_files(t_sample_name,
                                             n_sample_name,
                                             variant_call_tar,
                                             temp_dir)
        # write metrics to file
        metrics_file = write_qc_metric_to_file(
            t_bas_content,
            t_sample_name,
            n_bas_content,
            n_sample_name,
            genome_size,
            extracted_gender_file,
            extracted_purity_file,
            extracted_contamination_files,
            extracted_snv_file,
            extracted_indel_file,
            extracted_sv_file,
            extracted_cnv_file,
            temp_dir)
        # remove files that are not required in the output
        # clean_temp_dir(temp_dir, t_sample_name, n_sample_name)
        # tar all files in temp_dir to the ourput_tar
        try:
            with tarfile.open(output_tar, 'w:gz') as tar:
                tar.add(metrics_file, arcname=os.path.basename(metrics_file))
                for a_file in extracted_genotype_files:
                    tar.add(a_file, arcname=os.path.basename(a_file))
        except Exception as exc:
            raise RuntimeError('failed to create the final output: %s' % str(exc))
    logger.info('completed')


def extract_and_place_required_files(
        tumour_sample_name,
        normal_sample_name,
        variant_call_tar,
        temp_dir):
    # get a list of files that are required
    gender_file, purity_file, contamination_files, genotyping_files \
        = get_metrics_file_names(tumour_sample_name, normal_sample_name)

    if count_variants:
        snv_file, indel_file, sv_file, cnv_file \
            = get_variant_file_names(tumour_sample_name, normal_sample_name)

    # extract files from tar tar ball
    with tarfile.open(variant_call_tar, 'r:gz') as tar:
        logger.info('getting file list info from tar file %s', variant_call_tar)
        all_files = tar.getmembers()
        required_list = [gender_file] + [purity_file] + contamination_files + genotyping_files
        if count_variants:
            required_list.extend([snv_file, indel_file, sv_file, cnv_file])
        required_files = [a_file for a_file in all_files if a_file.name in required_list]
        # Check if all required files are found
        if len(required_files) < len(required_list):
            found_files = [a_file.name for a_file in required_files]
            missed_files = [a_file for a_file in required_list if a_file not in found_files]
            raise RuntimeError(
                'required files are not found in the variant call result tar ball: %s'
                % ', '.join(missed_files))
        logger.info('extracting files')
        tar.extractall(path=temp_dir, members=required_files)
        logger.info('extration done!')

    # move genotyping files to the top level within the temp_dir
    for a_file in genotyping_files:
        os.rename(
            os.path.join(temp_dir, a_file),
            os.path.join(temp_dir, os.path.basename(a_file)),
        )

    if count_variants:
        extracted_variant_files = [
            os.path.join(temp_dir, snv_file),
            os.path.join(temp_dir, indel_file),
            os.path.join(temp_dir, sv_file),
            os.path.join(temp_dir, cnv_file)]
    else:
        extracted_variant_files = [None, None, None, None]

    # return contamination files and gender check file
    return (
        os.path.join(temp_dir, gender_file),
        os.path.join(temp_dir, purity_file),
        [os.path.join(temp_dir, a_file) for a_file in contamination_files],
        [os.path.join(temp_dir, os.path.basename(a_file)) for a_file in genotyping_files],
        *extracted_variant_files
    )


def write_qc_metric_to_file(t_bas_content,
                            tumour_sample_name,
                            n_bas_content,
                            normal_sample_name,
                            genome_size,
                            extracted_gender_file,
                            extracted_purity_file,
                            extracted_contamination_files: list,
                            extracted_snv_file,
                            extracted_indel_file,
                            extracted_sv_file,
                            extracted_cnv_file,
                            temp_dir):

    t_con_file, n_con_file = extracted_contamination_files
    gender, f_m_gender, f_m_geno = get_gender_info_from_file(
        tumour_sample_name, extracted_gender_file)

    # try to created the output file
    output_file = os.path.join(
        temp_dir,
        f'{tumour_sample_name}_vs_{normal_sample_name}.ppcg_sanger_metrics.txt')
    if os.path.exists(output_file):
        raise RuntimeError('output file: %s already exists.' % output_file)

    # extract info from tumour_bas
    tumour_qc_metrics = [
        tumour_sample_name,
        'Tumour',
        get_seq_depth_from_bas(t_bas_content, genome_size),
        get_mapping_rate_from_bas(t_bas_content),
        get_average_insert_size_from_bas(t_bas_content),
        get_average_insert_size_sd_from_bas(t_bas_content),
        get_gc_r1_from_bas(t_bas_content),
        get_gc_r2_from_bas(t_bas_content),
        get_duplicate_r_rate_from_bas(t_bas_content),
        get_mismatched_pair_rate_from_bas(t_bas_content),
        get_contamination_from_file(t_con_file, tumour_sample_name),
        gender, f_m_gender, f_m_geno,
        get_purity_from_file(extracted_purity_file)
        ]

    # extract info from normal_bas
    normal_qc_metrics = [
        normal_sample_name,
        'Normal',
        get_seq_depth_from_bas(n_bas_content, genome_size),
        get_mapping_rate_from_bas(n_bas_content),
        get_average_insert_size_from_bas(n_bas_content),
        get_average_insert_size_sd_from_bas(n_bas_content),
        get_gc_r1_from_bas(n_bas_content),
        get_gc_r2_from_bas(n_bas_content),
        get_duplicate_r_rate_from_bas(n_bas_content),
        get_mismatched_pair_rate_from_bas(n_bas_content),
        get_contamination_from_file(n_con_file, normal_sample_name),
        'NA',
        'NA',
        'NA',
        'NA'
        ]

    output_header = OUTPUT_HEADER[:]

    if count_variants:
        # count variants
        tumour_qc_metrics.extend([get_snv_count(extracted_snv_file),
                                 get_indel_count(extracted_indel_file),
                                 get_sv_count(extracted_sv_file),
                                 get_cnv_count(extracted_cnv_file)])
        normal_qc_metrics.extend(['NA', 'NA', 'NA', 'NA'])
        output_header.extend(VARIANT_COUNT_HEADER)

    # write to the output
    try:
        with open(output_file, 'w') as out:
            out.write('\t'.join(output_header) + '\n')
            out.write('\t'.join([str(ele) for ele in tumour_qc_metrics]) + '\n')
            out.write('\t'.join([str(ele) for ele in normal_qc_metrics]) + '\n')
    except Exception as exc:
        raise RuntimeError('failed to write to output file: %s' % str(exc))
    return output_file


def get_abs_path(path):
    return_v = path
    if not os.path.isabs(path):
        return_v = os.path.abspath(path)
    return return_v


def check_file_exists(file_path):
    if not os.path.exists(file_path):
        raise RuntimeError('file %s does not exist.' % file_path)


def validate_bas(bas_file):
    file_name = os.path.basename(bas_file)
    if re.match(r'.+\.bam\.bas$', file_name) is None:
        raise RuntimeError('invalid BAS filename: %s, expecting ".bam.bas" suffix.' % bas_file)
    total_lines = 0
    with open(bas_file, 'r') as f:
        lines = f.readlines()
    if len(lines) <= 1:
        raise RuntimeError('too few lines in %s.' % bas_file)
    bas = [line.rstrip('\n').split('\t') for line in lines]
    if bas[0] != BAS_HEADER:
        raise RuntimeError('invalid BAS header in %s.' % bas_file)
    col_number = len(bas[0])
    for row in bas:
        if len(row) != col_number:
            raise RuntimeError('invalid row in BAS file: %s.' % '\t'.join(row))


def validate_tar_name(tar_file_name):
    if re.match(r'.+\.tar\.gz$', os.path.basename(tar_file_name)) is None:
        raise RuntimeError('%s should have ".tar.gz" suffix' % tar_file_name)


def get_bas_content(bas_file):
    with open(bas_file, 'r') as f:
        lines = f.readlines()
    bas = [line.rstrip('\n').split('\t') for line in lines]
    bas.pop(0)
    return bas


def get_sample_name_from_bas(bas_content):
    sample_names = list(set([rg[BAS_HEADER.index('sample')] for rg in bas_content]))
    if len(sample_names) > 1:
        raise RuntimeError(
            'invalid BAS file: too many sample names - sample names found %s' %
            ', '.join(sample_names))
    else:
        return sample_names[0]


def get_seq_depth_from_bas(bas_content, genome_size):
    total_mapped = sum([int(rg[BAS_HEADER.index('#_mapped_bases')]) for rg in bas_content])
    return total_mapped/genome_size


def get_mapping_rate_from_bas(bas_content):
    mapped_reads = sum([int(rg[BAS_HEADER.index('#_mapped_reads')]) for rg in bas_content])
    total_reads = sum([int(rg[BAS_HEADER.index('#_total_reads')]) for rg in bas_content])
    return mapped_reads/total_reads


def get_average_insert_size_from_bas(bas_content):
    insert_sizes = [float(rg[BAS_HEADER.index('mean_insert_size')]) for rg in bas_content]
    return sum(insert_sizes)/len(insert_sizes)


def get_average_insert_size_sd_from_bas(bas_content):
    insert_sizes = [float(rg[BAS_HEADER.index('insert_size_sd')]) for rg in bas_content]
    return sum(insert_sizes)/len(insert_sizes)


def get_gc_r1_from_bas(bas_content):
    total_r1_gc_b = sum([int(rg[BAS_HEADER.index('#_gc_bases_r1')]) for rg in bas_content])
    total_r1_b = sum(
        [int(rg[BAS_HEADER.index('#_total_reads_r1')]) * int(rg[BAS_HEADER.index('read_length_r1')])
            for rg in bas_content]
    )
    return total_r1_gc_b/total_r1_b


def get_gc_r2_from_bas(bas_content):
    total_r2_gc_b = sum([int(rg[BAS_HEADER.index('#_gc_bases_r2')]) for rg in bas_content])
    total_r2_b = sum(
        [int(rg[BAS_HEADER.index('#_total_reads_r2')]) * int(rg[BAS_HEADER.index('read_length_r2')])
            for rg in bas_content]
    )
    return total_r2_gc_b/total_r2_b


def get_duplicate_r_rate_from_bas(bas_content):
    total_duplicate_reads = sum(
        [int(rg[BAS_HEADER.index('#_duplicate_reads')]) for rg in bas_content])
    total_reads = sum([int(rg[BAS_HEADER.index('#_total_reads')]) for rg in bas_content])
    return total_duplicate_reads/total_reads


def get_mismatched_pair_rate_from_bas(bas_content):
    total_mapped_pairs = sum([int(rg[BAS_HEADER.index('#_mapped_pairs')]) for rg in bas_content])
    total_pairs = sum([int(rg[BAS_HEADER.index('#_total_reads_r1')]) for rg in bas_content])
    return 1-(total_mapped_pairs/total_pairs)


def get_contamination_from_file(con_file, sample_name):
    try:
        con_dict = json.loads(open(con_file).read())
    except Exception as exc:
        raise RuntimeError('can not load to dict: %s' % str(exc))
    try:
        con = con_dict[sample_name]['contamination']
    except Exception as exc:
        raise RuntimeError('can not find contamination value: %s' % str(exc))
    return con


def get_gender_info_from_file(tumour_sample_name, extraced_gender_file):
    try:
        gen_dict = json.loads(open(extraced_gender_file).read())
    except Exception as exc:
        raise RuntimeError('can not load to dict: %s' % str(exc))
    try:
        gender = gen_dict['tumours'][0]['gender']['gender']
        f_matched = gen_dict['tumours'][0]['gender']['frac_match_gender']
        f_matched_geno = gen_dict['tumours'][0]['genotype']['frac_matched_genotype']
        f_sample_name = gen_dict['tumours'][0]['sample']
    except Exception as exc:
        raise RuntimeError('can not find gender info: %s' % str(exc))
    if tumour_sample_name != f_sample_name:
        raise RuntimeError(
            'sample name does not match exptected in gender info file, expected %s, found %s'
            % (tumour_sample_name, f_sample_name))
    return (gender, f_matched, f_matched_geno)


def get_purity_from_file(purity_file):
    row_name = normal_contamination = None
    try:
        for line in open(purity_file, 'r').readlines():
            if re.match(r'^NormalContamination\s', line):
                row_name, normal_contamination = line.rstrip('\n').split(' ')
                break
    except Exception as exc:
        raise RuntimeError('can not load purity_file: %s' % str(exc))
    if not row_name:
        raise RuntimeError(
            'purity_file: %s does not have a line starting with "NormalContamination ".'
            % purity_file)
    if normal_contamination == '':
        raise RuntimeError(
            'normal contamination value in purity_file is invalid, expect a float, found nothing.')
    return normal_contamination


def get_metrics_file_names(tumour_sample_name, normal_sample_name):
    gender_file = f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/genotyped/result.json'
    contamination_files = [
        f'WGS_{tumour_sample_name}/contamination/result.json',
        f'WGS_{normal_sample_name}/contamination/result.json'
    ]
    purity_file = \
        f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/ascat/{tumour_sample_name}.samplestatistics.txt'
    genotyping_files = [
        f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/genotyped/{tumour_sample_name}.full_gender.tsv',
        f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/genotyped/{tumour_sample_name}.full_genotype.tsv',
        f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/genotyped/{normal_sample_name}.full_gender.tsv',
        f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/genotyped/{normal_sample_name}.full_genotype.tsv'
    ]
    return (
        gender_file,
        purity_file,
        contamination_files,
        genotyping_files
    )


def get_variant_file_names(tumour_sample_name, normal_sample_name):
    snv_file = f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/caveman/{tumour_sample_name}_vs_{normal_sample_name}.flagged.muts.vcf.gz'
    indel_file = f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/pindel/{tumour_sample_name}_vs_{normal_sample_name}.flagged.vcf.gz'
    sv_file = f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/brass/{tumour_sample_name}_vs_{normal_sample_name}.annot.vcf.gz'
    cnv_file = f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/ascat/{tumour_sample_name}.copynumber.caveman.vcf.gz'
    return (snv_file, indel_file, sv_file, cnv_file)


def get_snv_count(extracted_snv_file):
    check_file_exists(extracted_snv_file)
    cmd_return, return_code = exec_subp_and_wait(
        f'gunzip -c {extracted_snv_file} | grep -c "\tPASS\t"')
    # when no pattern is found, grep exit 1 but it should be a expected behavior
    if return_code == 1 and cmd_return != '0':
        raise RuntimeError(f"Subprocess failed with return code {returncode}!")
    return cmd_return


def get_indel_count(extracted_indel_file):
    check_file_exists(extracted_indel_file)
    cmd_return, return_code = exec_subp_and_wait(
        f'gunzip -c {extracted_indel_file} | grep -c "\tPASS\t"')
    # when no pattern is found, grep exit 1 but it should be a expected behavior
    if return_code == 1 and cmd_return != '0':
        raise RuntimeError(f"Subprocess failed with return code {returncode}!")
    return cmd_return


def get_sv_count(extracted_sv_file):
    check_file_exists(extracted_sv_file)
    cmd_return, return_code = exec_subp_and_wait(
        f'gunzip -c {extracted_sv_file}  | grep -v \'#\' | grep -c \'BAS=\'')
    # when no pattern is found, grep exit 1 but it should be a expected behavior
    if return_code == 1 and cmd_return != '0':
        raise RuntimeError(f"Subprocess failed with return code {returncode}!")
    if int(cmd_return) % 2 != 0:
        raise RuntimeError('counted BRASS calls returns an odd number.')
    # integer devision
    return str(int(cmd_return)//2)


def get_cnv_count(extracted_cnv_file):
    check_file_exists(extracted_cnv_file)
    count = 0
    with gzip.open(extracted_cnv_file, 'rt') as file:
        for line in file:
            if re.match(r'^#', line):
                continue
            else:
                columns = line.rstrip('\n').split('\t')
                if columns[9] != columns[10]:
                    count += 1
    return str(count)


def clean_temp_dir(temp_dir, tumour_sample_name, normal_sample_name):
    try:
        shutil.rmtree(os.path.join(temp_dir, f'WGS_{tumour_sample_name}_vs_{normal_sample_name}'))
        shutil.rmtree(os.path.join(temp_dir, f'WGS_{tumour_sample_name}'))
        shutil.rmtree(os.path.join(temp_dir, f'WGS_{normal_sample_name}'))
    except Exception as exc:
        raise RuntimeError('failed to remove temp_dir: %s' % str(exc))
