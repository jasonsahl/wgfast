#!/usr/bin/env python

"""test each function in the wg-fast
code"""

import unittest
import os
import tempfile
import shutil
import re
import sys
from wg_fast.util import *

curr_dir=os.getcwd()

class Test1(unittest.TestCase):
    def test_get_seq_name_basic_function(self):
        self.assertEqual(get_seq_name("/path/to/test.fasta"), "test.fasta")
    """tests the condition where you use a tilda instead of full path"""
    def test_get_seq_name_tilda(self):
        self.assertEqual(get_seq_name("~/test.fasta"), "test.fasta")
    """tests the case where no path is passed"""
    def test_get_seq_name_empty(self):
        self.assertEqual(get_seq_name(""), "")
    """tests the case where something weird is passed"""
    def test_get_seq_name_wrong_slash(self):
        self.assertEqual(get_seq_name("\wrong\way"), "\\wrong\\way")

class Test2(unittest.TestCase):
    def test_get_readFile_components_basic_function(self):
        self.assertEqual(get_readFile_components("/path/to/file.gz"), ('/path/to', 'file', '.gz'))
    def test_get_readFile_components_tilda(self):
        self.assertEqual(get_readFile_components("~/path/to/file.gz"), ('~/path/to', 'file', '.gz'))
    def test_get_readFile_components_non_gz(self):
        self.assertEqual(get_readFile_components("~/path/to/file.fasta"), ('~/path/to', 'file', '.fasta'))
    def test_get_readFile_components_wrong_slash(self):
        self.assertEqual(get_readFile_components("~\path\to\file.fasta"), ('', '~\\path\to\x0cile', '.fasta'))
    def test_get_readFile_components_empty(self):
        self.assertEqual(get_readFile_components(""), ('', '', ''))

class Test5(unittest.TestCase):
    def test_sort_information_basic_function(self):
        self.assertEqual(sort_information("ADK1::460"), 460)
    def test_sort_information_cant_parse(self):
        self.assertRaises(TypeError, sort_information, "ADK1__460")
    def test_sort_information_no_input(self):
        self.assertRaises(TypeError, sort_information, None)

class Test6(unittest.TestCase):
    def test_matrix_to_fasta_basic_function(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered.vcf")
        fp = open(fpath, "w")
        fp.write("LocusID\tReference\tgenome1\tgenome2\n")
        fp.write("ADK::1\tA\tT\tT\n")
        fp.write("ADK::2\tT\tT\tT\n")
        fp.close()
        self.assertEqual(matrix_to_fasta(fpath, "%s/outfile.txt" % tdir), [">ReferenceAT", ">genome1TT", ">genome2TT"])
        shutil.rmtree(tdir)
    def test_matrix_to_fasta_unequal_fields(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered.vcf")
        fp = open(fpath, "w")
        fp.write("LocusID\tReference\tgenome1\tgenome2\n")
        fp.write("ADK::1\tA\tT\tT\n")
        fp.write("ADK::2\tT\tT\n")
        fp.close()
        self.assertEqual(matrix_to_fasta(fpath, "%s/outfile.txt" % tdir), [">ReferenceAT", ">genome1TT"])
        shutil.rmtree(tdir)
    def test_matrix_to_fasta_multiple_states(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered.vcf")
        fp = open(fpath, "w")
        fp.write("LocusID\tReference\tgenome1\tgenome2\n")
        fp.write("ADK::1\tA\tT\tT\n")
        fp.write("ADK::2\tT\tT\tTT\n")
        fp.close()
        self.assertEqual(matrix_to_fasta(fpath, "%s/outfile.txt" % tdir), [">ReferenceAT", ">genome1TT", ">genome2TTT"])
        shutil.rmtree(tdir)

class Test7(unittest.TestCase):
    def test_write_reduced_matrix_basic_function(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered.vcf")
        fp = open(fpath, "w")
        fp.write("LocusID\tReference\tgenome1\tgenome2\t#SNPcall\n")
        fp.write("ADK::1\tA\tT\tT\n")
        fp.write("ADK::2\tT\tT\tT\n")
        fp.close()
        self.assertEqual(write_reduced_matrix(fpath), [4, 4])
        shutil.rmtree(tdir)
    def test_write_reduced_matrix_odd_field_numbers(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered.vcf")
        fp = open(fpath, "w")
        fp.write("LocusID\tReference\tgenome1\tgenome2\t#SNPcall\n")
        fp.write("ADK::1\tA\tT\tT\n")
        fp.write("ADK::2\tT\tT\n")
        fp.close()
        self.assertEqual(write_reduced_matrix(fpath), [4, 3])
        shutil.rmtree(tdir)

class Test8(unittest.TestCase):
    def test_make_temp_matrix_basic_function(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.matrix")
        fp = open(fpath, "w")
        fp.write("LocusID\tReference\tgenome1\tgenome2\n")
        fp.write("ADK::1\tA\tT\tT\n")
        fp.write("ADK::2\tT\tT\tT\n")
        fp.write("ADK::3\tG\tG\tT\n")
        fp.close()
        vpath = os.path.join(tdir,"testfile.filtered.vcf")
        vp = open(vpath, "w")
        vp.write("ADK::1\tG\n")
        vp.close()
        self.assertEqual(make_temp_matrix(vpath,fpath,"test"), {'ADK::1': 'G', 'ADK::2': 'N', 'ADK::3': 'N'})
        shutil.rmtree(tdir)
    def test_make_temp_matrix_no_matches(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.matrix")
        fp = open(fpath, "w")
        fp.write("LocusID\tReference\tgenome1\tgenome2\n")
        fp.write("ADK::1\tA\tT\tT\n")
        fp.write("ADK::2\tT\tT\tT\n")
        fp.write("ADK::3\tG\tG\tT\n")
        fp.close()
        vpath = os.path.join(tdir,"testfile.filtered.vcf")
        vp = open(vpath, "w")
        vp.write("ADK::6\tA\n")
        vp.close()
        self.assertEqual(make_temp_matrix(vpath,fpath,"test"), {'ADK::1': 'N', 'ADK::2': 'N', 'ADK::3': 'N', 'ADK::6': 'A'})
        shutil.rmtree(tdir)

class Test9(unittest.TestCase):
    def test_grab_names_basic_function(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"test.filtered.vcf")
        fpath2 = os.path.join(tdir,"name_1.filtered.vcf")
        os.chdir("%s" % tdir)
        fp = open(fpath, "w")
        fp2 = open(fpath2, "w")
        fp.close()
        fp2.close()
        self.assertEqual(grab_names(), ['name_1', 'test'])
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)

#class Test10(unittest.TestCase):
#    def test_process_coverage_basic_function(self):
#        tdir = tempfile.mkdtemp(prefix="filetest_",)
#        fpath = os.path.join(tdir,"ECOLI_coverage.sample_summary")
#        os.chdir("%s" % tdir)
#        fp = open(fpath,"w")
#        fp.write("sample_id       total   mean    granular_third_quartile granular_median granular_first_quartile %_bases_above_15\n")
#        fp.write("ECOLI    2050    3.82    6       5       4       0.0\n")
#        fp.write("Total    2050    3.82    N/A     N/A     N/A")
#        fp.close()
#        self.assertEqual(process_coverage("ECOLI"), {'ECOLI':'3.82'})
#        os.chdir("%s" % curr_dir)
#        shutil.rmtree(tdir)
#    def test_process_coverage_missing_match(self):
#        tdir = tempfile.mkdtemp(prefix="filetest_",)
#        fpath = os.path.join(tdir,"ECOLI_coverage.sample_summary")
#        os.chdir("%s" % tdir)
#        fp = open(fpath,"w")
#        fp.write("sample_id       total   mean    granular_third_quartile granular_median granular_first_quartile %_bases_above_15\n")
#        fp.write("EOLI    2050    3.82    6       5       4       0.0\n")
#        fp.write("Total    2050    3.82    N/A     N/A     N/A")
#        fp.close()
#        self.assertRaises(TypeError, process_coverage, "ECOLI")
#        os.chdir("%s" % curr_dir)
#        shutil.rmtree(tdir)

class Test11(unittest.TestCase):
    def test_find_two_report_error(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"all_patristic_distances.txt")
        os.chdir("%s" % tdir)
        fp = open(fpath,"w")
        fp.write("Distance between 'E2348_69_allexternalnucmer' and 'H10407_allexternalnucmer': 1.39030683167\n")
        fp.close()
        self.assertRaises(TypeError, find_two_new, ['E2348_69_allexternalnucmer'], ['H10407_allexternalnucmer','E2348_69_allexternalnucmer'])
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)
    def test_find_two_basic_function(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"all_patristic_distances.txt")
        os.chdir("%s" % tdir)
        fp = open(fpath,"w")
        fp.write("Distance between 'E2348_69_allexternalnucmer' and 'H10407_allexternalnucmer': 1.39030683167\n")
        fp.write("Distance between 'E2348_69_allexternalnucmer' and 'O157_H7_sakai_allexternalnucmer': 4.53192608862\n")
        fp.write("Distance between 'E2348_69_allexternalnucmer' and 'Reference': 1.29611949657")
        fp.close()
        self.assertEqual(find_two_new(fpath, ["E2348_69_allexternalnucmer"]),((('E2348_69_allexternalnucmer', 'O157_H7_sakai_allexternalnucmer', '4.53192608862'), ('E2348_69_allexternalnucmer', 'H10407_allexternalnucmer', '1.39030683167')),(('Reference', 'E2348_69_allexternalnucmer', '1.29611949657'),)))
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)
    def test_find_two_reversed(self):
        tdir = tempfile.mkdtemp(prefix="reversed_",)
        fpath = os.path.join(tdir,"all_patristic_distances_reversed.txt")
        os.chdir("%s" % tdir)
        fp = open(fpath,"w")
        fp.write("Distance between 'E2348_69_allexternalnucmer' and 'H10407_allexternalnucmer': 1.39030683167\n")
        fp.write("Distance between 'E2348_69_allexternalnucmer' and 'O157_H7_sakai_allexternalnucmer': 4.53192608862\n")
        fp.write("Distance between 'Reference' and 'E2348_69_allexternalnucmer': 0.0941892612547")
        fp.close()
        self.assertEqual(find_two_new(fpath, ["E2348_69_allexternalnucmer"]),((('E2348_69_allexternalnucmer', 'O157_H7_sakai_allexternalnucmer', '4.53192608862'), ('E2348_69_allexternalnucmer', 'H10407_allexternalnucmer', '1.39030683167')),(('Reference', 'E2348_69_allexternalnucmer', '0.0941892612547'),)))
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)
class Test12(unittest.TestCase):
    def test_get_closest_dists_basic_function(self):
        self.assertEqual(get_closest_dists_new((('ECOLI_ISO2', 'SSON_046_allexternalnucmer', '0.08198920048'), ('ECOLI_ISO2', 'H10407_allexternalnucmer', '1.3087194675e-06')),['ECOLI', 'ECOLI_IS03_L007', 'ECOLI_ISO2']),(['SSON_046_allexternalnucmer0.08198920048', 'H10407_allexternalnucmer1.3087194675e-06']))

class Test13(unittest.TestCase):
    def test_calculate_pairwise_tree_dists_basic_function(self):
        """distances were taken directly from Dendropy, run outside of the pipeline"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"test.tree")
        os.chdir("%s" % tdir)
        fp = open(fpath,"w")
        fp.write("((H10407_all:0.00000095154411648831,SSON_046_all:0.08004748973386473232):0.09609983129754622044,(E2348_69_all:1.63492542157114217893,O157_H7_sakai_all:4.52711175943011490119):0.00000095154411648831,Reference:0.00000095154411648831):0.0;")
        fp.close()
        self.assertEqual(calculate_pairwise_tree_dists(fpath, "tmp.out"),[0.08004844127798122, 1.7310271559569212, 4.623213493815894, 0.0961017343857792, 1.8110736941466694, 4.703260032005642, 0.17614827257552745, 6.162037181001257, 1.634927324659375, 4.527113662518348])
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)

class Test14(unittest.TestCase):
    def test_subsample_snps(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"test.matrix")
        os.chdir("%s" % tdir)
        fp = open(fpath,"w")
        fp.write("LocusID\tReference\tgenome1\tgenome2\n")
        fp.write("ADK::1\tA\tT\tT\n")
        fp.write("ADK::2\tT\tT\tT\n")
        fp.write("ADK::3\tG\tG\tT\n")
        fp.close()
        self.assertEqual(subsample_snps(fpath,{'ECOLI': ['SSON_046_allexternalnucmer', 'H10407_allexternalnucmer', 'E2348_69_allexternalnucmer', 'O157_H7_sakai_allexternalnucmer'], 'ECOLI_ISO2': ['H10407_allexternalnucmer', 'SSON_046_allexternalnucmer', 'E2348_69_allexternalnucmer', 'O157_H7_sakai_allexternalnucmer'], 'ECOLI_IS03_L007': ['SSON_046_allexternalnucmer', 'H10407_allexternalnucmer', 'E2348_69_allexternalnucmer', 'O157_H7_sakai_allexternalnucmer']}, {'ECOLI': 11, 'ECOLI_ISO2': 14, 'ECOLI_IS03_L007': 11},4),['ADK::1','ADK::2','ADK::3'])
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)

class Test15(unittest.TestCase):
    def test_process_temp_matrices_dev(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir, "tmp.tree")
        fpath2 = os.path.join(tdir, "all.distances")
        os.chdir("%s" % tdir)
        fp = open(fpath,"w")
        fp2 = open(fpath2, "w")
        fp.write("(E2348_69_allexternalnucmer:1.29611826169204280568,O157_H7_sakai_allexternalnucmer:3.23580782692392832089,(Reference:0.00000061743893499046,((QUERY___ECOLI_ISO2:0.000001,H10407_allexternalnucmer:0.00000030871946749523):0.00000030871946749523,((QUERY___ECOLI:0.000001,QUERY___ECOLI_IS03_L007:0.000001):0.0,SSON_046_allexternalnucmer:0.04099394588025907088):0.04099394588025907088):0.09418733509629113876):0.00000061743893499046);")
        fp.close()
        fp2.write("Distance between 'E2348_69_allexternalnucmer' and 'O157_H7_sakai_allexternalnucmer': 4.53192608862\n")
        fp2.write("Distance between 'E2348_69_allexternalnucmer' and 'Reference': 1.29611949657\n")
        fp2.write("Distance between 'E2348_69_allexternalnucmer' and 'QUERY___ECOLI_ISO2': 1.39030752295\n")
        fp2.write("Distance between 'E2348_69_allexternalnucmer' and 'H10407_allexternalnucmer': 1.39030683167\n")
        fp2.write("Distance between 'E2348_69_allexternalnucmer' and 'QUERY___ECOLI': 1.43130116011\n")
        fp2.write("Distance between 'E2348_69_allexternalnucmer' and 'QUERY___ECOLI_IS03_L007': 1.43130116011\n")
        fp2.write("Distance between 'E2348_69_allexternalnucmer' and 'SSON_046_allexternalnucmer': 1.47229410599\n")
        fp2.write("Distance between 'O157_H7_sakai_allexternalnucmer' and 'Reference': 3.2358090618\n")
        fp2.write("Distance between 'O157_H7_sakai_allexternalnucmer' and 'QUERY___ECOLI_ISO2': 3.32999708818\n")
        fp2.write("Distance between 'O157_H7_sakai_allexternalnucmer' and 'H10407_allexternalnucmer': 3.3299963969\n")
        fp2.write("Distance between 'O157_H7_sakai_allexternalnucmer' and 'QUERY___ECOLI': 3.37099072534\n")
        fp2.write("Distance between 'O157_H7_sakai_allexternalnucmer' and 'QUERY___ECOLI_IS03_L007': 3.37099072534\n")
        fp2.write("Distance between 'O157_H7_sakai_allexternalnucmer' and 'SSON_046_allexternalnucmer': 3.41198367122\n")
        fp2.write("Distance between 'Reference' and 'QUERY___ECOLI_ISO2': 0.0941892612547\n")
        fp2.write("Distance between 'Reference' and 'H10407_allexternalnucmer': 0.0941885699742\n")
        fp2.write("Distance between 'Reference' and 'QUERY___ECOLI': 0.135182898415\n")
        fp2.write("Distance between 'Reference' and 'QUERY___ECOLI_IS03_L007': 0.135182898415\n")
        fp2.write("Distance between 'Reference' and 'SSON_046_allexternalnucmer': 0.176175844296\n")
        fp2.write("Distance between 'QUERY___ECOLI_ISO2' and 'H10407_allexternalnucmer': 1.3087194675e-06\n")
        fp2.write("Distance between 'QUERY___ECOLI_ISO2' and 'QUERY___ECOLI': 0.0409962545997\n")
        fp2.write("Distance between 'QUERY___ECOLI_ISO2' and 'QUERY___ECOLI_IS03_L007': 0.0409962545997\n")
        fp2.write("Distance between 'QUERY___ECOLI_ISO2' and 'SSON_046_allexternalnucmer': 0.08198920048\n")
        fp2.write("Distance between 'H10407_allexternalnucmer' and 'QUERY___ECOLI': 0.0409955633192\n")
        fp2.write("Distance between 'H10407_allexternalnucmer' and 'QUERY___ECOLI_IS03_L007': 0.0409955633192\n")
        fp2.write("Distance between 'H10407_allexternalnucmer' and 'SSON_046_allexternalnucmer': 0.0819885091995\n")
        fp2.write("Distance between 'QUERY___ECOLI' and 'QUERY___ECOLI_IS03_L007': 2e-06\n")
        fp2.write("Distance between 'QUERY___ECOLI' and 'SSON_046_allexternalnucmer': 0.0409949458803\n")
        fp2.write("Distance between 'QUERY___ECOLI_IS03_L007' and 'SSON_046_allexternalnucmer': 0.0409949458803")
        fp2.close()
        #self.assertEqual(process_temp_matrices(({'ECOLI': ['SSON_046_allexternalnucmer', 'H10407_allexternalnucmer', 'E2348_69_allexternalnucmer', 'O157_H7_sakai_allexternalnucmer'], 'ECOLI_ISO2': ['H10407_allexternalnucmer', 'SSON_046_allexternalnucmer', 'E2348_69_allexternalnucmer', 'O157_H7_sakai_allexternalnucmer'], 'ECOLI_IS03_L007': ['SSON_046_allexternalnucmer', 'H10407_allexternalnucmer', 'E2348_69_allexternalnucmer', 'O157_H7_sakai_allexternalnucmer']},fpath,2,fpath2),(

class Test16(unittest.TestCase):
    def test_get_all_snps(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered.vcf")
        fp = open(fpath, "w")
        fp.write("LocusID\tReference\tgenome1\tgenome2\n")
        fp.write("ADK::1\tA\tT\tT\n")
        fp.write("ADK::2\tT\tT\n")
        fp.close()
        self.assertEqual(get_all_snps(fpath), ['ADK::1','ADK::2'])
        shutil.rmtree(tdir)

class Test17(unittest.TestCase):
    def test_find_used_snps_basic(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered.vcf")
        fp = open(fpath, "w")
        fp.write("gi|16120353|ref|NC_003143.1|::35501\tG\n")
        fp.write("gi|16120353|ref|NC_003143.1|::52924\tC\n")
        fp.write("gi|16120353|ref|NC_003143.1|::55551\tN\n")
        fp.close()
        os.chdir("%s" % tdir)
        self.assertEqual(find_used_snps(),{'testfile':2})
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)
    def test_find_used_snps_all_missing(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered.vcf")
        fp = open(fpath, "w")
        fp.write("gi|16120353|ref|NC_003143.1|::35501\tN\n")
        fp.write("gi|16120353|ref|NC_003143.1|::52924\tN\n")
        fp.write("gi|16120353|ref|NC_003143.1|::55551\tN\n")
        fp.close()
        os.chdir("%s" % tdir)
        self.assertEqual(find_used_snps(),{'testfile':0})
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)
    def test_find_used_snps_blank_line(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered.vcf")
        fp = open(fpath, "w")
        fp.write("gi|16120353|ref|NC_003143.1|::35501\tG\n")
        fp.write("gi|16120353|ref|NC_003143.1|::52924\tC\n")
        fp.write("gi|16120353|ref|NC_003143.1|::55551\t-\n")
        fp.write(" n")
        fp.close()
        os.chdir("%s" % tdir)
        self.assertRaises(TypeError, find_used_snps, ())
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)

class Test18(unittest.TestCase):
    def test_branch_lengths_to_decimals(self):
        self.assertEqual(branch_lengths_2_decimals("('E2348_69_allexternalnucmer':1.63492542157,'O157_H7_sakai_allexternalnucmer':4.52711175943):9.51544116488e-07,Reference:9.51544116488e-07)"),("('E2348_69_allexternalnucmer':1.634925,'O157_H7_sakai_allexternalnucmer':4.527112):0.000001,Reference:0.000001);"))

class Test19(unittest.TestCase):
    def test_parse_likelihoods(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.filtered.vcf")
        fp = open(fpath, "w")
        fp.write("ECOLI\tI6\t0.298716\t0.298716\n")
        fp.write("ECOLI\tI5\t0.298714\t0.597430\n")
        fp.write("ECOLI\tI4\t0.298714\t0.896143\n")
        fp.close()
        self.assertEqual(parse_likelihoods(fpath), {'ECOLI':['0.298716','0.298714','0.298714']})
        shutil.rmtree(tdir)

class Test20(unittest.TestCase):
    def test_fasta_to_tab(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.fasta")
        fp = open(fpath, "w")
        fp.write(">id\n")
        fp.write("ATGC")
        fp.close()
        self.assertEqual(fasta_to_tab(fpath), ["id","ATGC"])
        shutil.rmtree(tdir)

class Test21(unittest.TestCase):
    def test_tab_to_fasta(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.tab")
        fp = open(fpath, "w")
        fp.write(">id\tATCG\n")
        fp.close()
        self.assertEqual(tab_to_fasta(fpath), [">id","ATCG"])
        shutil.rmtree(tdir)

class Test22(unittest.TestCase):
    def test_tab_to_matrix(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.tab")
        fp = open(fpath, "w")
        fp.write("id\tATCG\n")
        fp.write("id2\tTTCG\n")
        fp.close()
        self.assertEqual(tab_to_matrix(fpath), [['id', 'id2'], ['A', 'T'], ['T', 'T'], ['C', 'C'], ['G', 'G']])
        shutil.rmtree(tdir)

class Test23(unittest.TestCase):
    def test_prune_fasta(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.tab")
        fp = open(fpath, "w")
        fp.write(">id\nATCG\n")
        fp.write(">id2\nTTCG\n")
        fp.close()
        self.assertEqual(prune_fasta("id",fpath,"%s/outfile.txt" % tdir), ['id2'])
        shutil.rmtree(tdir)

class Test24(unittest.TestCase):
    def test_filter_alignment(self):
        """tests to make sure that this function
        throws out all characters that it is supposed to"""
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        fpath = os.path.join(tdir,"testfile.tab")
        fp = open(fpath, "w")
        fp.write("throw away line\n")
        fp.write("X\tA\tN\tX\t-\n")
        fp.close()
        self.assertEqual(filter_alignment(fpath), [['A']])
        shutil.rmtree(tdir)

class Test25(unittest.TestCase):
    def test_compare_subsample_results(self):
        tdir = tempfile.mkdtemp(prefix="filetest_",)
        os.chdir("%s" % tdir)
        fpath = os.path.join(tdir,"sample..testfile.tab")
        fp = open(fpath, "w")
        fp.write("Cluster0	Cluster0	100.00	15	0	0	1	15	1	15	1e-07	30.2")
        fp.close()
        os.chdir("%s" % curr_dir)
        shutil.rmtree(tdir)
        os.system("rm tab_matrix tab.filtered out.tab out.fasta")

if __name__ == "__main__":
    unittest.main()
    main()
