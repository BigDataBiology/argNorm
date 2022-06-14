import unittest
import argNorm.normalizers as argnorm

class TestARGSOAPNormalizer(unittest.TestCase):

    def test_raw_reads_mode_input(self):
        self.__test_run(is_hamronized=False, mode='reads')

    def test_raw_orfs_mode_input(self):
        self.__test_run(is_hamronized=False, mode='orfs')

    def test_ham_reads_mode_input(self):
        self.__test_run(is_hamronized=True, mode='reads')
        
    def test_ham_orfs_mode_input(self):
        self.__test_run(is_hamronized=True, mode='orfs')

    def __test_run(self, is_hamronized, mode):
        if is_hamronized:
            prefix = 'examples/hamronized/'
        else:
            prefix = 'examples/raw/'
        try:
            norm = argnorm.ARGSOAPNormalizer(is_hamronized=is_hamronized, mode=mode)
            norm.run(input_file=prefix + f'args-oap.sarg.{mode}.tsv')
            raised = False
        except:
            raised = True
        self.assertEqual(raised, False)