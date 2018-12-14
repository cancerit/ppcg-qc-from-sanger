import os
import shutil
import sys
import errno
import tarfile
import logging
import json
import re
import gzip
from statistics import median
from tempfile import TemporaryDirectory
from . import exec_subp_and_wait, check_file_exists, format_float
from typing import Dict, List, Any


# setup logs
logger = logging.getLogger('sanger_qc_extractor')
# create console handler and set level to debug
cha = logging.StreamHandler()
cha.setLevel(logging.DEBUG)
cha.setFormatter(logging.Formatter('%(asctime)s %(levelname)8s - %(message)s',
                                   '%Y-%m-%d %H:%M:%S'))
# add ch to logger
logger.addHandler(cha)
logger.setLevel(logging.DEBUG)


def set_extractor_logger_level(logging_level):
    global logger
    logger.setLevel(logging_level)


def join_and_median(a_list):
    return ','.join(
        [format_float(a_number) for a_number in a_list]
    ), format_float(median(a_list))


class SangerQcMetricsExtractor(object):

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
        '#_inter_chr_pairs'
    ]

    VALID_METADATA_KEYS = [
        'donor_id',
        'donor_uuid',
        'tumour_sequencing_year',
        'tumour_sequencer',
        'tumour_id',
        'tumour_uuid',
        'normal_sequencing_year',
        'normal_sequencer',
        'normal_id',
        'normal_uuid'
    ]

    NA_STRING = 'NA'

    def __init__(self, tumour_bas, normal_bas, genome_size, variant_call_tar, output_dir, count_variants, metadata: Dict[str, str] = None):
        '''
        the main function for handling the whole QC metrics extraction process
        '''
        if not isinstance(genome_size, int):
            raise RuntimeError('genome_size is not int')

        # if bas files hava valid columns
        self.validate_bas(tumour_bas)
        self.validate_bas(normal_bas)
        self.validate_tar_name(variant_call_tar)

        self.genome_size = genome_size
        self.variant_call_tar = variant_call_tar
        self.t_bas_content = self.get_bas_content(tumour_bas)
        self.t_sample_name = self.get_sample_name_from_bas(self.t_bas_content)
        self.n_bas_content = self.get_bas_content(normal_bas)
        self.n_sample_name = self.get_sample_name_from_bas(self.n_bas_content)
        self.t_rg_ids = self.get_rg_ids_from_bas(self.t_bas_content)
        self.n_rg_ids = self.get_rg_ids_from_bas(self.n_bas_content)
        self.count_variants = count_variants
        self.output_dir = output_dir

        # extracted required files and move genotyping files to output folder
        self.extracted_files: Dict[str, Any] = self.extract_and_place_required_files()

        self.metadata = {}
        if isinstance(metadata, dict):
            self.metadata = self.validate_metadata(metadata)

    def extract_and_place_required_files(self):
        # get a list of files that are required
        gender_file, purity_file, contamination_files, genotyping_files = \
            self.get_metrics_file_names(self.t_sample_name, self.n_sample_name)

        if self.count_variants:
            snv_file, indel_file, sv_file, cnv_file = \
               self.get_variant_file_names(self.t_sample_name, self.n_sample_name)

        # extract files from tar tar ball
        with tarfile.open(self.variant_call_tar, 'r:gz') as tar:
            logger.info('getting file list info from tar file %s', self.variant_call_tar)
            all_files = tar.getmembers()
            required_list = [gender_file] + [purity_file] + contamination_files + genotyping_files
            if self.count_variants:
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
            tar.extractall(path=self.output_dir, members=required_files)
            logger.info('extration done!')

        # move genotyping files to the top level within the output_dir
        for a_file in genotyping_files:
            os.rename(
                os.path.join(self.output_dir, a_file),
                os.path.join(self.output_dir, os.path.basename(a_file)),
            )

        if self.count_variants:
            extracted_variant_files = {
                'snv': os.path.join(self.output_dir, snv_file),
                'indel': os.path.join(self.output_dir, indel_file),
                'sv': os.path.join(self.output_dir, sv_file),
                'cnv': os.path.join(self.output_dir, cnv_file)
            }
        else:
            extracted_variant_files = {
                'snv': None,
                'indel': None,
                'sv': None,
                'cnv': None
            }

        # return contamination files and gender check file
        return {
            'gender': os.path.join(self.output_dir, gender_file),
            'purity': os.path.join(self.output_dir, purity_file),
            'contamination': [os.path.join(self.output_dir, a_file) for a_file in contamination_files],
            'genotyping': [os.path.join(self.output_dir, os.path.basename(a_file)) for a_file in genotyping_files],
            'variants': extracted_variant_files
        }

    def get_metrics(self) -> List[str]:
        gender, f_m_gender, f_m_geno = self.get_gender_info_from_file()

        to_return = [
            self.t_sample_name,
            self.metadata.get('tumour_id', SangerQcMetricsExtractor.NA_STRING),
            self.metadata.get('tumour_uuid', SangerQcMetricsExtractor.NA_STRING),
            self.metadata.get('tumour_sequencing_year', SangerQcMetricsExtractor.NA_STRING),
            self.metadata.get('tumour_sequencer', SangerQcMetricsExtractor.NA_STRING),
            ','.join(self.t_rg_ids),
            *self.get_tumour_seq_depth(),
            *self.get_tumour_mapping_rate(),
            *self.get_tumour_insert_sizes(),
            *self.get_tumour_insert_size_sds(),
            *self.get_tumour_gc_r1(),
            *self.get_tumour_gc_r2(),
            *self.get_tumour_duplicate_r_rate(),
            *self.get_tumour_mismatched_pair_rate(),
            *self.get_tumour_contamination(),
            gender, f_m_gender, f_m_geno,
            self.get_purity_from_file(),
            self.n_sample_name,
            self.metadata.get('normal_id', SangerQcMetricsExtractor.NA_STRING),
            self.metadata.get('normal_uuid', SangerQcMetricsExtractor.NA_STRING),
            self.metadata.get('normal_sequencing_year', SangerQcMetricsExtractor.NA_STRING),
            self.metadata.get('normal_sequencer', SangerQcMetricsExtractor.NA_STRING),
            ','.join(self.n_rg_ids),
            *self.get_normal_seq_depth(),
            *self.get_normal_mapping_rate(),
            *self.get_normal_insert_sizes(),
            *self.get_normal_insert_size_sds(),
            *self.get_normal_gc_r1(),
            *self.get_normal_gc_r2(),
            *self.get_normal_duplicate_r_rate(),
            *self.get_normal_mismatched_pair_rate(),
            *self.get_normal_contamination(),
            self.metadata.get('donor_id', SangerQcMetricsExtractor.NA_STRING),
            self.metadata.get('donor_uuid', SangerQcMetricsExtractor.NA_STRING)
        ]

        if self.count_variants:
            # count variants
            to_return.extend([
                '\t'.join(self.get_snv_count()),
                '\t'.join(self.get_indel_count()),
                '\t'.join(self.get_sv_count()),
                self.get_cnv_count()
            ])

        return to_return

    def get_genotyping_files(self):
        return self.extracted_files['genotyping']

    def clean_output_dir(self):
        try:
            shutil.rmtree(os.path.join(self.output_dir, f'WGS_{self.t_sample_name}_vs_{self.n_sample_name}'))
            shutil.rmtree(os.path.join(self.output_dir, f'WGS_{self.t_sample_name}'))
            shutil.rmtree(os.path.join(self.output_dir, f'WGS_{self.n_sample_name}'))
        except Exception as exc:
            raise RuntimeError('failed to remove temp_dir: %s' % str(exc))

    def get_gender_info_from_file(self):
        try:
            gen_dict = json.loads(open(self.extracted_files['gender']).read())
        except Exception as exc:
            raise RuntimeError('can not load to dict: %s' % str(exc))
        try:
            gender = gen_dict['tumours'][0]['gender']['gender']
            f_matched = gen_dict['tumours'][0]['gender']['frac_match_gender']
            f_matched_geno = gen_dict['tumours'][0]['genotype']['frac_matched_genotype']
            f_sample_name = gen_dict['tumours'][0]['sample']
        except Exception as exc:
            raise RuntimeError('can not find gender info: %s' % str(exc))
        if self.t_sample_name != f_sample_name:
            raise RuntimeError(
                'sample name does not match exptected in gender info file, expected %s, found %s'
                % (self.t_sample_name, f_sample_name))
        return (gender, f_matched, f_matched_geno)

    def get_purity_from_file(self):
        purity_file = self.extracted_files['purity']
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

    def get_snv_count(self):
        return self.get_v_count(self.extracted_files['variants']['snv'])

    def get_indel_count(self):
        return self.get_v_count(self.extracted_files['variants']['indel'])

    def get_sv_count(self):
        extracted_sv_file = self.extracted_files['variants']['sv']
        check_file_exists(extracted_sv_file)

        sv_per_chr = {}
        with gzip.open(extracted_sv_file, 'rt') as v:
            for line in v:
                if not line.startswith('#'):
                    ele = line.rstrip('\n').split('\t')
                    if not ele[0] in sv_per_chr:
                        sv_per_chr[ele[0]] = [0, 0]
                    sv_per_chr[ele[0]][1] += 1
                    if re.search(r'BAS=', ele[7]):
                        sv_per_chr[ele[0]][0] += 1

        count_filtered = sum([n[0] for n in sv_per_chr.values()])
        count_all = sum([n[1] for n in sv_per_chr.values()])

        if int(count_all) % 2 != 0 or int(count_filtered) % 2 != 0:
            raise RuntimeError('counted BRASS calls returns an odd number.')
        # integer devision 
        return f'{str(int(count_filtered)//2)}/{str(int(count_all)//2)}', json.dumps({k:[n/2 for n in v] for k,v in sv_per_chr.items()})

    def get_cnv_count(self):
        extracted_cnv_file = self.extracted_files['variants']['cnv']
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

    @staticmethod
    def validate_bas(bas_file):

        check_file_exists(bas_file)

        file_name = os.path.basename(bas_file)
        if re.match(r'.+\.bam\.bas$', file_name) is None:
            raise RuntimeError('invalid BAS filename: %s, expecting ".bam.bas" suffix.' % bas_file)
        with open(bas_file, 'r') as f:
            lines = f.readlines()
        if len(lines) <= 1:
            raise RuntimeError('too few lines in %s.' % bas_file)
        bas = [line.rstrip('\n').split('\t') for line in lines]
        if bas[0] != SangerQcMetricsExtractor.BAS_HEADER:
            raise RuntimeError('invalid BAS header in %s.' % bas_file)
        col_number = len(bas[0])
        for row in bas:
            if len(row) != col_number:
                raise RuntimeError('invalid row in BAS file: %s.' % '\t'.join(row))

    @staticmethod
    def validate_tar_name(tar_file_name):
        if re.match(r'.+\.tar\.gz$', os.path.basename(tar_file_name)) is None:
            raise RuntimeError('%s should have ".tar.gz" suffix' % tar_file_name)

    @staticmethod
    def get_bas_content(bas_file):
        with open(bas_file, 'r') as f:
            lines = f.readlines()
        bas = [line.rstrip('\n').split('\t') for line in lines]
        bas.pop(0)
        return bas

    @staticmethod
    def get_rg_ids_from_bas(bas_content):
        return [
                line[SangerQcMetricsExtractor.BAS_HEADER.index('readgroup')] for line in bas_content
            ]

    @staticmethod
    def get_sample_name_from_bas(bas_content):
        sample_names = list(set([rg[SangerQcMetricsExtractor.BAS_HEADER.index('sample')] for rg in bas_content]))
        if len(sample_names) > 1:
            raise RuntimeError(
                'invalid BAS file: too many sample names - sample names found %s' %
                ', '.join(sample_names))
        else:
            return sample_names[0]

    @staticmethod
    def get_seq_depth_from_bas(bas_content, genome_size):
        mapped_per_lane = [int(rg[SangerQcMetricsExtractor.BAS_HEADER.index('#_mapped_bases')]) for rg in bas_content]
        total_mapped = sum(mapped_per_lane)
        return ','.join([format_float(mapped/genome_size) for mapped in mapped_per_lane]), format_float(total_mapped/genome_size)

    def get_tumour_seq_depth(self):
        return self.get_seq_depth_from_bas(self.t_bas_content, self.genome_size)

    def get_normal_seq_depth(self):
        return self.get_seq_depth_from_bas(self.n_bas_content, self.genome_size)

    @staticmethod
    def get_mapping_rate_from_bas(bas_content):
        mapped_reads = [int(rg[SangerQcMetricsExtractor.BAS_HEADER.index('#_mapped_reads')]) for rg in bas_content]
        total_reads = [int(rg[SangerQcMetricsExtractor.BAS_HEADER.index('#_total_reads')]) for rg in bas_content]
        return join_and_median(
            [m_read/t_read for m_read, t_read in zip(mapped_reads, total_reads)]
        )

    def get_tumour_mapping_rate(self):
        return self.get_mapping_rate_from_bas(self.t_bas_content)

    def get_normal_mapping_rate(self):
        return self.get_mapping_rate_from_bas(self.n_bas_content)

    @staticmethod
    def get_insert_sizes_from_bas(bas_content):
        insert_sizes = [float(rg[SangerQcMetricsExtractor.BAS_HEADER.index('mean_insert_size')]) for rg in bas_content]
        return join_and_median(insert_sizes)

    def get_tumour_insert_sizes(self):
        return self.get_insert_sizes_from_bas(self.t_bas_content)

    def get_normal_insert_sizes(self):
        return self.get_insert_sizes_from_bas(self.n_bas_content)

    @staticmethod
    def get_insert_size_sds_from_bas(bas_content):
        insert_size_sds = [float(rg[SangerQcMetricsExtractor.BAS_HEADER.index('insert_size_sd')]) for rg in bas_content]
        return join_and_median(insert_size_sds)

    def get_tumour_insert_size_sds(self):
        return self.get_insert_size_sds_from_bas(self.t_bas_content)

    def get_normal_insert_size_sds(self):
        return self.get_insert_size_sds_from_bas(self.n_bas_content)

    @staticmethod
    def get_gc_from_bas(bas_content, r1=True):
        '''
        will return GC content of mate 1 in read pairs by defult, if r1=False, return GC content of mate 2 in read pairs
        '''

        bas_gc_index = SangerQcMetricsExtractor.BAS_HEADER.index('#_gc_bases_r1')
        bas_tb_index = SangerQcMetricsExtractor.BAS_HEADER.index('#_total_reads_r1')
        bas_len_index = SangerQcMetricsExtractor.BAS_HEADER.index('read_length_r1')

        if not r1:
            bas_gc_index = SangerQcMetricsExtractor.BAS_HEADER.index('#_gc_bases_r2')
            bas_tb_index = SangerQcMetricsExtractor.BAS_HEADER.index('#_total_reads_r2')
            bas_len_index = SangerQcMetricsExtractor.BAS_HEADER.index('read_length_r2')

        r_gc_bs = [int(rg[bas_gc_index]) for rg in bas_content]
        r_bs = [
                int(rg[bas_tb_index]) * int(rg[bas_len_index])
                for rg in bas_content
            ]

        return join_and_median(
            [r_gc_b/r_b for r_gc_b, r_b in zip(r_gc_bs, r_bs)]
        )

    def get_tumour_gc_r1(self):
        return self.get_gc_from_bas(self.t_bas_content)

    def get_normal_gc_r1(self):
        return self.get_gc_from_bas(self.n_bas_content)

    def get_tumour_gc_r2(self):
        return self.get_gc_from_bas(self.t_bas_content, r1=False)

    def get_normal_gc_r2(self):
        return self.get_gc_from_bas(self.n_bas_content, r1=False)

    @staticmethod
    def get_duplicate_r_rate_from_bas(bas_content):
        duplicate_reads = [int(rg[SangerQcMetricsExtractor.BAS_HEADER.index('#_duplicate_reads')]) for rg in bas_content]
        reads = [int(rg[SangerQcMetricsExtractor.BAS_HEADER.index('#_total_reads')]) for rg in bas_content]

        return join_and_median(
            [d_r/r for d_r, r in zip(duplicate_reads, reads)]
        )

    def get_tumour_duplicate_r_rate(self):
        return self.get_duplicate_r_rate_from_bas(self.t_bas_content)

    def get_normal_duplicate_r_rate(self):
        return self.get_duplicate_r_rate_from_bas(self.n_bas_content)

    @staticmethod
    def get_mismatched_pair_rate_from_bas(bas_content):
        mapped_pairs = [int(rg[SangerQcMetricsExtractor.BAS_HEADER.index('#_mapped_pairs')]) for rg in bas_content]
        pairs = [int(rg[SangerQcMetricsExtractor.BAS_HEADER.index('#_total_reads_r1')]) for rg in bas_content]

        return join_and_median(
            [1-(m_p/p) for m_p, p in zip(mapped_pairs, pairs)]
        )

    def get_tumour_mismatched_pair_rate(self):
        return self.get_mismatched_pair_rate_from_bas(self.t_bas_content)

    def get_normal_mismatched_pair_rate(self):
        return self.get_mismatched_pair_rate_from_bas(self.n_bas_content)

    @staticmethod
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

    @staticmethod
    def get_variant_file_names(tumour_sample_name, normal_sample_name):
        snv_file = f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/caveman/{tumour_sample_name}_vs_{normal_sample_name}.flagged.muts.vcf.gz'
        indel_file = f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/pindel/{tumour_sample_name}_vs_{normal_sample_name}.flagged.vcf.gz'
        sv_file = f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/brass/{tumour_sample_name}_vs_{normal_sample_name}.annot.vcf.gz'
        cnv_file = f'WGS_{tumour_sample_name}_vs_{normal_sample_name}/ascat/{tumour_sample_name}.copynumber.caveman.vcf.gz'
        return (snv_file, indel_file, sv_file, cnv_file)

    @staticmethod
    def get_v_count(v_file):
        count_per_chr={}
        with gzip.open(v_file, 'rt') as v:
            for line in v:
                if not line.startswith('#'):
                    ele = line.rstrip('\n').split('\t')
                    if ele[0] not in count_per_chr:
                        count_per_chr[ele[0]] = [0, 0]

                    count_per_chr[ele[0]][1] += 1
                    if ele[6] == 'PASS':
                        count_per_chr[ele[0]][0] += 1
        count_filtered = sum([c_t[0] for c_t in count_per_chr.values()])
        count_all = sum([c_t[1] for c_t in count_per_chr.values()])

        return f'{count_filtered}/{count_all}', json.dumps(count_per_chr)

    @staticmethod
    def get_contamination_from_file(rg_ids, con_file, sample_name):
        try:
            con_dict = json.loads(open(con_file).read())
        except Exception as exc:
            raise RuntimeError('can not load to dict: %s' % str(exc))
        try:
            contamination_per_lane = []
            for rg_id in rg_ids:
                contamination_per_lane.append(con_dict[sample_name]['by_readgroup'][rg_id]['contamination'])
        except Exception as exc:
            raise RuntimeError('can not find contamination value: %s' % str(exc))

        # the contaminations can be really small (eg: 0.00000001), thus not a good idea to use format_float
        return ','.join([str(cpl) for cpl in contamination_per_lane]), median(contamination_per_lane)

    def get_tumour_contamination(self):
        return self.get_contamination_from_file(self.t_rg_ids, self.extracted_files['contamination'][0], self.t_sample_name)

    def get_normal_contamination(self):
        return self.get_contamination_from_file(self.n_rg_ids, self.extracted_files['contamination'][1], self.n_sample_name)

    @staticmethod
    def validate_metadata(metadata: Dict[str, str]):
        valid_metadata = {}
        available = sorted(
            list(
                set(SangerQcMetricsExtractor.VALID_METADATA_KEYS) & set(metadata.keys())
            ),
            key=lambda x: SangerQcMetricsExtractor.VALID_METADATA_KEYS.index(x)
        )
        if available:
            logger.debug('%s are found in the given metadata', ', '.join(available))

        for key in available:
            if isinstance(metadata[key], str):
                valid_metadata[key] = metadata[key]
            else:
                logger.warning('%s is in metadata but not a string, skip.', key)

        return valid_metadata
