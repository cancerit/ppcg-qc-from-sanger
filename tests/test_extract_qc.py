import os
import unittest
import pytest
from uuid import uuid4
from tempfile import TemporaryDirectory
from argparse import Namespace
import tarfile
import ppcg_qc_from_sanger.extract_qc as extract_qc

class PpcgQcFromSangerExtractQcTestCase(unittest.TestCase):

    def setUp(self):
        self.test_dir = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'data')

    def test_get_purity(self):
        testing_function = extract_qc.get_purity_from_file
        with TemporaryDirectory() as temp_dir:
            test_file_content = 'NormalContamination 0.002\n'
            test_file = self.write_to_file(temp_dir, test_file_content)
            purity = testing_function(test_file)
            self.assertEqual(purity, '0.002')

            with pytest.raises(RuntimeError) as exc:
                purity = testing_function('some_file')
            self.assertIn('can not load purity_file', str(exc.value))

            test_file_content = 'NormalContamination \n'
            test_file = self.write_to_file(temp_dir, test_file_content)
            with pytest.raises(RuntimeError) as exc:
                purity = testing_function(test_file)
            self.assertIn('first line in purity_file is invalid', str(exc.value))

            test_file_content = 'normalContamination \n'
            test_file = self.write_to_file(temp_dir, test_file_content)
            with pytest.raises(RuntimeError) as exc:
                purity = testing_function(test_file)
            self.assertIn('first line in purity_file is invalid', str(exc.value))

    def write_to_file(self, dir, content, file_name=None):
        if file_name is None:
            test_file = os.path.join(dir, str(uuid4()))
        else:
            test_file = os.path.join(dir, file_name)
        with open(test_file, 'w') as f:
            f.write(content)
        return test_file

    def test_validate_tar_name(self):
        testing_function = extract_qc.validate_tar_name
        with pytest.raises(RuntimeError) as exc:
            result = testing_function(os.path.join('x', 'y', 'z', '.tar.gz'))
        self.assertIn('should have ".tar.gz" suffix', str(exc.value))
        with pytest.raises(RuntimeError) as exc:
            result = testing_function(os.path.join('x', 'y', 'z.gz'))
        self.assertIn('should have ".tar.gz" suffix', str(exc.value))

    def test_validate_bas(self):
        testing_function = extract_qc.validate_bas

        # test invalid file name
        with pytest.raises(RuntimeError) as exc:
            test = testing_function('test.bam.basss')
        self.assertIn('invalid BAS filename', str(exc.value))

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

    def test_extract_and_place_required_files(self):
        testing_function = extract_qc.extract_and_place_required_files
        with TemporaryDirectory() as temp_dir:
            e_gender_file, e_purity_file, e_contamination_files, e_genotype_files = \
                extract_qc.get_metrics_file_names('tumour_sample', 'normal_sample')
            return_files = testing_function(
                'tumour_sample',
                'normal_sample',
                os.path.join(self.test_dir, 'test.tar.gz'),
                temp_dir)
            self.assertEqual(
                (
                    os.path.join(temp_dir, e_gender_file),
                    os.path.join(temp_dir, e_purity_file),
                    [os.path.join(temp_dir, a_file) for a_file in e_contamination_files],
                    [os.path.join(temp_dir, os.path.basename(a_file)) for a_file in e_genotype_files]
                ),
                return_files
            )
            for a_file in e_genotype_files:
                self.assertTrue(os.path.exists(
                    os.path.join(temp_dir, os.path.basename(a_file)))
                )

            # test wrong sample name
            with pytest.raises(RuntimeError) as exc:
                return_files = testing_function(
                    'tumour_sample1',
                    'normal_sample',
                    os.path.join(self.test_dir, 'test.tar.gz'),
                    temp_dir)
            self.assertIn('required files are not found', str(exc.value))

    def test_extract_from_sanger(self):
        '''
        testing the whole process
        '''
        testing_function = extract_qc.extract_from_sanger
        with TemporaryDirectory() as temp_dir:
            good_args = Namespace(
                tumour_bas=os.path.join(self.test_dir, 'tumour_sample.bam.bas'),
                normal_bas=os.path.join(self.test_dir, 'normal_sample.bam.bas'),
                variant_call_tar=os.path.join(self.test_dir, 'test.tar.gz'),
                output_tar=os.path.join(temp_dir, 'result.tar.gz'),
                genome_size=3137454505,
                debug=False
            )
            self.assertEqual(testing_function(good_args), None)
            # check if the metrics are as expected.
            metrics_file = 'tumour_sample_vs_normal_sample.ppcg_sanger_metrics.txt'
            with tarfile.open(good_args.output_tar, 'r:gz') as tar:
                result_files = [a_file.name for a_file in tar.getmembers()]
                print(result_files)
                tar.extract(metrics_file, path=temp_dir)
            self.assertEqual(
                open(os.path.join(temp_dir, metrics_file), 'r').read(),
                open(os.path.join(self.test_dir, 'expected_metrics.txt')).read()
            )
            # check if expected files are in the tar
            e_genotype_files = \
                extract_qc.get_metrics_file_names('tumour_sample', 'normal_sample')[3]
            for a_file in e_genotype_files:
                print(a_file)
                print(result_files)
                self.assertTrue(os.path.basename(a_file) in result_files)

            # test wrong genome size type
            # a shallow copy of good_args
            args = Namespace(**vars(good_args))
            args.genome_size = '1'
            with pytest.raises(RuntimeError) as exc:
                testing_function(args)
            self.assertEqual('genome_size is not int', str(exc.value))

            # test non exist file
            args = Namespace(**vars(good_args))
            args.tumour_bas = 'some should not exist file.file'
            with pytest.raises(RuntimeError) as exc:
                testing_function(args)
            self.assertIn('does not exist', str(exc.value))
