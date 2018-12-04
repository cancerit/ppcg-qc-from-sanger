import os
import sys
import tarfile
import logging
import re
import json
from tempfile import TemporaryDirectory
from .sanger_qc_extractor import SangerQcMetricsExtractor, set_extractor_logger_level
from . import check_file_exists_for_user
from typing import List, Dict, Tuple, Any

# TODO Change inputs to lists, complain if tumour samples have duplicated names, complain if missing bas files, warning if there're more tumour or normal bas files!
# TODO use tumour sample names as folder names in output, put genotype files in it
# TODO all metrics to one file
# TODO tar the folder and metrics file in one ball!!


# setup logs
logger = logging.getLogger('ppcg_qc_from_sanger')
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
    'Tumour sample name',
    'Tumour sample ID',
    'Tumour sample UUID',
    'Tumour sequencing year',
    'Tumour sequencer',
    'Tumour ReadGroup IDs',
    'Tumour depth per RG',
    'Tumour total depth',
    'Tumour fraction of mapped reads per RG',
    'Tumour mean fraction of mapped reads',
    'Tumour insert size per RG',
    'Tumour mean Insert size',
    'Tumour insert size sd per RG',
    'Tumour mean insert size sd',
    'Tumour r1 GC content per RG',
    'Tumour mean r1 GC content',
    'Tumour r2 GC content per RG',
    'Tumour mean r2 GC content',
    'Tumour fraction of duplicated reads per RG',
    'Tumour mean fraction of duplicated reads',
    'Tumour fraction of mis-matched pairs per RG',
    'Tumour mean fraction of mis-matched pairs',
    'Tumour contamination per RG',
    'Tumour mean contamination',
    'Tumour sex',
    'Tumour fraction of matched sex with Normal',
    'Tumour fraction of matched genotype with Normal',
    'Normal contamination in Tumour',
    'Normal sample name',
    'Normal sample ID',
    'Normal sample UUID',
    'Normal sequencing year',
    'Normal sequencer',
    'Normal ReadGroup IDs',
    'Normal depth per RG',
    'Normal total depth',
    'Normal fraction of mapped reads per RG',
    'Normal mean fraction of mapped reads',
    'Normal insert size per RG',
    'Normal mean Insert size',
    'Normal insert size sd per RG',
    'Normal mean insert size sd',
    'Normal r1 GC content per RG',
    'Normal mean r1 GC content',
    'Normal r2 GC content per RG',
    'Normal mean r2 GC content',
    'Normal fraction of duplicated reads per RG',
    'Normal mean fraction of duplicated reads',
    'Normal fraction of mis-matched pairs per RG',
    'Normal mean fraction of mis-matched pairs',
    'Normal contamination per RG',
    'Normal mean contamination',
    'Donor ID',
    'Donor UUID']

VARIANT_COUNT_HEADER = [
    'Number of SNVs (PASS/All)', 'Number of INDELs (PASS/All)', 'Number of SVs (PASS/All)', 'Number of CNVs']

PPCG_META_TO_EXTRACTOR_MAP = {
    'donor_id': 'donor_id',
    'donor_uuid': 'donor_uuid',
    'sequencer': ['tumour_sequencing_year', 'normal_sequencing_year'],
    'sequencing_year': ['tumour_sequencing_year', 'normal_sequencing_year'],
    'sample_id': ['tumour_id', 'normal_id'],
    'sample_uuid': ['tumour_uuid', 'normal_uuid']
}


def extract_from_sanger(args):
    '''
    the main function for handling the whole QC metrics extraction process
    '''
    if not args.debug:
        logger.setLevel(logging.INFO)
        set_extractor_logger_level(logging.INFO)

    check_paras(args)

    genome_size = args.genome_size
    output_tar = os.path.abspath(args.output_tar)

    logger.info('checking tumour BAS file(s)..')
    tumour_bas = get_all_bas(args.tumour_bas)
    logger.info('checking normal BAS file(s)..')
    normal_bas = get_all_bas(args.normal_bas)

    if args.metadata:
        logger.info('checking metadata file(s)..')
        sample_id_meta, sample_uuid_meta = get_sample_meta(get_all_meta(args.metadata))

    logger.info('checking Sanger variant call result file(s)..')
    variant_call_tars = get_all_variant_call_tar(args.variant_call_tar)

    t_name_bas, n_name_bas, t_n_pair_tar = \
        get_validated_t_n_pair_and_bas_lists(tumour_bas, normal_bas, variant_call_tars)

    t_n_pair_meta = {}
    if args.metadata:
        logger.info('read in metadata file(s)..')
        t_n_pair_meta = get_t_n_pair_meta(sample_id_meta, sample_uuid_meta, t_n_pair_tar.keys())

    count_variants = args.count_variants
    # print('count_variants', count_variants)

    # create a temp dir to store extracted files
    with TemporaryDirectory() as temp_dir:
        output_metrics_file = os.path.join(temp_dir, 'ppcg_sanger_metrics.txt')
        genotyping_files = []

        with open(output_metrics_file, 'w') as o:
            header = OUTPUT_HEADER
            if count_variants:
                header += VARIANT_COUNT_HEADER
            o.write('\t'.join(header) + '\n')

            for t_n_pair, v_tar in t_n_pair_tar.items():
                logger.info('processing sanger call tar: %s', v_tar)
                logger.debug('processing sanger call tar: %s, tumour: %s, normal: %s, genome_size: %s, count_v: %s, meta_data: %s', v_tar, t_name_bas[t_n_pair[0]], n_name_bas[t_n_pair[1]], genome_size, count_variants, json.dumps(t_n_pair_meta.get(t_n_pair, {})))
                extractor = SangerQcMetricsExtractor(t_name_bas[t_n_pair[0]], n_name_bas[t_n_pair[1]], genome_size, v_tar, temp_dir, count_variants, t_n_pair_meta.get(t_n_pair, None))
                o.write('\t'.join(extractor.get_metrics()) + '\n')
                genotyping_files.extend(extractor.get_genotyping_files())
                extractor.clean_output_dir()

        # tar all files in temp_dir to the ourput_tar
        try:
            with tarfile.open(output_tar, 'w:gz') as tar:
                tar.add(output_metrics_file, arcname=os.path.basename(output_metrics_file))
                for a_file in genotyping_files:
                    tar.add(a_file, arcname=os.path.basename(a_file))
        except Exception as exc:
            logger.critical('failed to create the final output: %s', str(exc))
            sys.exit(1)
    logger.info('completed')


def check_paras(args):
    if not isinstance(args.genome_size, int):
        logger.critical('genome_size is not int')
        sys.exit(1)

    output_tar = os.path.abspath(args.output_tar)
    # if tar file has '.tar.gz' extension
    SangerQcMetricsExtractor.validate_tar_name(output_tar)

    # test if output is writable
    if not os.path.exists(output_tar):
        try:
            with open(output_tar, 'w') as out:
                out.write('place holder\n')
        except OSError as exc:
            logger.critical('output is not writable: %s.', str(exc))
            sys.exit(1)
        finally:
            os.remove(output_tar)
    else:
        logger.critical('existing output file: %s.', output_tar)
        sys.exit(1)


def append_to_file_path_list(a_path, path_list):
    if a_path in path_list:
        logger.warning(f'Duplicated input of a file: {a_path}, skip.')
    else:
        path_list.append(a_path)
    return path_list


def get_all_bas(input_abs: List[str]):
    to_return = []
    for path in input_abs:
        check_file_exists_for_user(path, logger)
        if os.path.isdir(path):
            logger.warning('%s is a directory, will take all BAS files in the folder, but not any file in a sub directory.', path)
            for a_file in os.scandir(path):
                if os.path.isfile(a_file) and re.match(r'.+\.bam\.bas$', a_file.name):
                    SangerQcMetricsExtractor.validate_bas(a_file)
                    to_return = append_to_file_path_list(os.path.abspath(a_file), to_return)
        else:
            SangerQcMetricsExtractor.validate_bas(path)
            to_return = append_to_file_path_list(os.path.abspath(path), to_return)
    return to_return


def get_all_variant_call_tar(input_call_tars: List[str]):
    to_return = []
    for path in input_call_tars:
        check_file_exists_for_user(path, logger)
        if os.path.isdir(path):
            logger.warning('%s is a directory, will take all tar.gz files in the folder, but not any file in a sub directory.', path)
            for a_file in os.scandir(path):
                if os.path.isfile(a_file) and re.match(r'.+\.tar\.gz$', a_file.name):
                    to_return = append_to_file_path_list(os.path.abspath(a_file), to_return)
        else:
            SangerQcMetricsExtractor.validate_tar_name(path)
            to_return = append_to_file_path_list(os.path.abspath(path), to_return)
    return to_return


def get_validated_t_n_pair_and_bas_lists(tumour_bas: List[str], normal_bas: List[str], variant_call_tars: List[str]) -> Tuple[Dict[str, str], Dict[str, str], Dict[Tuple[str, str], str]]:
    t_n_pair_tar: Dict[Tuple[str, str], str] = get_all_t_n_pairs(variant_call_tars)
    expected_tumours = [a_pair[0] for a_pair in t_n_pair_tar.keys()]
    expected_normals = [a_pair[1] for a_pair in t_n_pair_tar.keys()]
    t_name_bas: Dict[str, str] = get_sample_names_bas_file_dict(tumour_bas)
    n_name_bas: Dict[str, str] = get_sample_names_bas_file_dict(normal_bas)

    # if all expected tumour have bas
    not_found = sorted(set(expected_tumours) - set(t_name_bas.keys()))
    if not_found:
        logger.critical('Missing BAS files for tumour samples: %s', ', '.join(not_found))
        sys.exit(1)
    # if all expected normal have bas
    not_found = sorted(set(expected_normals) - set(n_name_bas.keys()))
    if not_found:
        logger.critical('Missing BAS files for normal samples: %s', ', '.join(not_found))
        sys.exit(1)

    return t_name_bas, n_name_bas, t_n_pair_tar


def get_all_t_n_pairs(variant_call_tars) -> Dict[Tuple[str, str], str]:
    t_n_pair_tar = {}
    for a_tar in variant_call_tars:
        t_name = n_name = None
        with tarfile.open(a_tar, 'r:gz') as tar:
            logger.info('validating tar file %s', a_tar)
            all_files = tar.getmembers()
            for a_file in all_files:
                matches = re.match(r'^WGS_([\w\-]+)_vs_([\w\-]+)$', a_file.name)
                if matches:
                    t_name = matches.group(1)
                    n_name = matches.group(2)
                    break
        if not t_name:
            logger.critical(f'Not a valid Sanger Variant Call result archive: {a_tar}')
            sys.exit(1)
        t_n_pair_tar[(t_name, n_name)] = a_tar
    return t_n_pair_tar


def get_sample_names_bas_file_dict(bas_list):
    return {
        SangerQcMetricsExtractor.get_sample_name_from_bas(SangerQcMetricsExtractor.get_bas_content(bas)): bas
        for bas in bas_list
    }


def get_all_meta(metadata_paths):
    to_return = []
    for path in metadata_paths:
        check_file_exists_for_user(path, logger)
        if os.path.isdir(path):
            logger.warning('%s is a directory, will take all tsv files in the folder, but not any file in a sub directory.', path)
            for a_file in os.scandir(path):
                if os.path.isfile(a_file) and re.match(r'.+\.tsv$', a_file.name):
                    to_return = append_to_file_path_list(os.path.abspath(a_file), to_return)
        else:
            to_return = append_to_file_path_list(os.path.abspath(path), to_return)
    return to_return


def get_t_n_pair_meta(sample_id_meta, sample_uuid_meta, t_n_pairs: List[str]) -> Tuple[Dict[str, dict], Dict[str, dict]]:

    to_return = {}

    for tumour, normal in t_n_pairs:
        # make sure it's a new copy
        t_dict = {**sample_uuid_meta.get(tumour, {})}
        if not t_dict:
            t_dict = {**sample_id_meta.get(tumour, {})}

        # make sure it's a new copy
        n_dict = {**sample_uuid_meta.get(normal, {})}
        if not n_dict:
            n_dict = {**sample_id_meta.get(normal, {})}

        tmp_dict = {}
        for a_key in ['donor_id', 'donor_uuid']:
            n_dict.pop(a_key, None)
            if a_key in t_dict:
                tmp_dict[PPCG_META_TO_EXTRACTOR_MAP[a_key]] = t_dict.pop(a_key)

        tmp_dict.update({
            PPCG_META_TO_EXTRACTOR_MAP[a_key][0]: t_dict.get(a_key)
            for a_key in t_dict
        })
        tmp_dict.update({
            PPCG_META_TO_EXTRACTOR_MAP[a_key][1]: n_dict.get(a_key)
            for a_key in n_dict
        })

        if tmp_dict:
            to_return[(tumour, normal)] = tmp_dict

    return to_return


def get_sample_meta(meta_files) -> Tuple[Dict[str, dict], Dict[str, dict]]:
    '''
    process tsv files into valid metadata
    '''
    sample_id_meta = {}
    sample_uuid_meta = {}

    for meta_file in meta_files:
        with open(meta_file, 'r') as meta:
            logger.debug('processing metadata file: %s...', meta_file)
            sample_id_index = sample_uuid_index = None
            lower_case_header = [ele.lower() for ele in meta.readline().split('\t')]
            try:
                sample_id_index = lower_case_header.index('sample_id')
            except ValueError:
                logger.warning('can not find Sample_ID column in input meta: %s', meta_file)

            try:
                sample_id_index = lower_case_header.index('sample_uuid')
            except ValueError:
                logger.warning('can not find Sample_UUID column in input meta: %s', meta_file)

            if not any([sample_id_index, sample_uuid_index]):
                logger.critical('can not find Sample_ID or Sample_UUID column in input meta: %s', meta_file)
                sys.exit(1)

            for line in meta:
                line_split = line.split('\t')
                valid_dict = get_valid_meta_dict(lower_case_header, line_split)
                if sample_id_index:
                    sample_id = line_split[sample_id_index]
                    logger.debug('found metadata for sample id: %s, valid meta: \n%s', sample_id, json.dumps(valid_dict))
                    if sample_id in sample_id_meta:
                        logger.critical('duplicated sample_id: %s is found in: %s', sample_id, meta_file)
                        sys.exit(1)
                    sample_id_meta[sample_id] = valid_dict
                if sample_uuid_index:
                    sample_uuid = line_split[sample_uuid_index]
                    logger.debug('found metadata for sample id: %s, valid meta: \n%s', sample_uuid, json.dumps(valid_dict))
                    if sample_uuid in sample_uuid_meta:
                        logger.critical('duplicated sample_uuid: %s is found in: %s', sample_uuid, meta_file)
                        sys.exit(1)
                    sample_uuid_meta[sample_uuid] = {**valid_dict}

    return sample_id_meta, sample_uuid_meta


def get_valid_meta_dict(header, line_split):
    to_return = {}
    for a_key in PPCG_META_TO_EXTRACTOR_MAP.keys():
        if a_key in header:
            to_return[a_key] = line_split[header.index(a_key)]
    return to_return
