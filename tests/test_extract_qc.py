import os
import unittest
import pytest
from uuid import uuid4
from tempfile import TemporaryDirectory
from argparse import Namespace
import tarfile
import ppcg_qc_from_sanger.extract_qc as extract_qc
from ppcg_qc_from_sanger.sanger_qc_extractor import SangerQcMetricsExtractor


class PpcgQcFromSangerExtractQcTestCase(unittest.TestCase):

    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')

    def write_to_file(self, dir, content, file_name=None):
        if file_name is None:
            test_file = os.path.join(dir, str(uuid4()))
        else:
            test_file = os.path.join(dir, file_name)
        with open(test_file, 'w') as f:
            f.write(content)
        return test_file

    def test_validate_tar_name(self):
        testing_function = SangerQcMetricsExtractor.validate_tar_name
        with pytest.raises(RuntimeError) as exc:
            result = testing_function(os.path.join('x', 'y', 'z', '.tar.gz'))
        self.assertIn('should have ".tar.gz" suffix', str(exc.value))
        with pytest.raises(RuntimeError) as exc:
            result = testing_function(os.path.join('x', 'y', 'z.gz'))
        self.assertIn('should have ".tar.gz" suffix', str(exc.value))

    def test_validate_bas(self):
        testing_function = SangerQcMetricsExtractor.validate_bas

        # test invalid file name
        with pytest.raises(RuntimeError) as exc:
            test = testing_function('test.bam.basss')
        self.assertIn('file test.bam.basss does not exist.', str(exc.value))

        # test valid content
        test_file = os.path.join(self.test_dir, 'test_sample.bam.bas')
        self.assertEqual(testing_function(test_file), None)

        test_content = open(test_file).read()
        with TemporaryDirectory() as temp_dir:

            # test bad header
            test_file_content = 'bas_file_name\n' + test_content
            test_file = self.write_to_file(temp_dir, test_file_content, 'test_1.bam.bas')
            with pytest.raises(RuntimeError) as exc:
                test = testing_function(test_file)
            self.assertIn('invalid BAS header', str(exc.value))

            # test fewer columns in a row
            test_file_content = test_content + '-\ttest_sample\n'
            test_file = self.write_to_file(temp_dir, test_file_content, 'test_2.bam.bas')
            with pytest.raises(RuntimeError) as exc:
                test = testing_function(test_file)
            self.assertIn('invalid row in BAS file', str(exc.value))

    def test_extract_from_sanger(self):
        '''
        testing the whole process
        '''
        testing_function = extract_qc.extract_from_sanger
        with TemporaryDirectory() as temp_dir:
            good_args = Namespace(
                tumour_bas=[os.path.join(self.test_dir, 'tumour_sample.bam.bas')],
                normal_bas=[os.path.join(self.test_dir, 'normal_sample.bam.bas')],
                variant_call_tar=[os.path.join(self.test_dir, 'test.tar.gz')],
                output_tar=os.path.join(temp_dir, 'result.tar.gz'),
                genome_size=3137454505,
                debug=False,
                count_variants=False,
                metadata=[]
            )
            # run the function to have outputs
            self.assertEqual(testing_function(good_args), None)

            # check if the metrics are as expected.
            metrics_file = 'ppcg_sanger_metrics.txt'
            with tarfile.open(good_args.output_tar, 'r:gz') as tar:
                result_files = [a_file.name for a_file in tar.getmembers()]
                tar.extract(metrics_file, path=temp_dir)
            self.assertEqual(
                open(os.path.join(temp_dir, metrics_file), 'r').read(),
                open(os.path.join(self.test_dir, 'expected_metrics.txt')).read()
            )

            # check if expected files are in the tar
            e_genotype_files = \
                SangerQcMetricsExtractor.get_metrics_file_names('tumour_sample', 'normal_sample')[3]
            for a_file in e_genotype_files:
                self.assertTrue(os.path.basename(a_file) in result_files)

            # test wrong genome size type
            # a shallow copy of good_args
            args = Namespace(**vars(good_args))
            args.genome_size = '1'
            with pytest.raises(SystemExit) as exc:
                testing_function(args)
            self.assertEqual(exc.type, SystemExit)

            # test non exist file
            args = Namespace(**vars(good_args))
            args.tumour_bas = 'some should not exist file.file'
            with pytest.raises(SystemExit) as exc:
                testing_function(args)
            self.assertEqual(exc.type, SystemExit)

            # test if variants can be counted correctly
            args = Namespace(**vars(good_args))
            args.count_variants = True
            args.output_tar = os.path.join(temp_dir, 'result_2.tar.gz')
            # run the function to have outputs
            self.assertEqual(testing_function(args), None)
            metrics_file = 'ppcg_sanger_metrics.txt'
            with tarfile.open(args.output_tar, 'r:gz') as tar:
                result_files = [a_file.name for a_file in tar.getmembers()]
                tar.extract(metrics_file, path=temp_dir)
            self.assertEqual(
                open(os.path.join(temp_dir, metrics_file), 'r').read(),
                open(os.path.join(self.test_dir, 'expected_metrics_with_variant_counts.txt')).read()
            )
