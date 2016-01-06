#!/usr/bin/python
# coding: utf-8

import unittest
import os

from lib import ApytramNeeds
from lib import BlastPlus
from lib import Trinity
from lib import Aligner

def search(program):
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    
    if fpath:
        if is_exe(program):	   
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


class TestConfigApytram(unittest.TestCase):

	#Check Blastp:
	def test_BlastP(self):
		#print "\nSearch for blastp"
		Search = search("blastp")
		#print Search
		self.assertIsNotNone(Search)

        #Check makeblastdb :
        def test_makeblastdb(self):
                #print "\nSearch for makeblastdb"
                Search = search("makeblastdb")
                #print Search
                self.assertIsNotNone(Search)

        #Check blastdbcmd:
        def test_blastdbcmd(self):
                #print "\nSearch for blastdbcmd"
                Search = search("blastdbcmd")
                #print Search
                self.assertIsNotNone(Search)

        #Check exonerate:
        def test_exonerate(self):
                #print "\nSearch for exonerate"
                Search = search("exonerate")
                #print Search
                self.assertIsNotNone(Search)

        #Check Trinity:
        def test_Trinity(self):
                #print "\nSearch for Trinity"
                Search = search("Trinity")
                #print Search
                self.assertIsNotNone(Search)
		
        #Check samtools:
        def test_samtools(self):
                #print "\nSearch for samtools"
                Search = search("samtools")
                #print Search
		#if Search:
		#	print "        !! Check samtools version = 0.1.19"
                self.assertIsNotNone(Search)

	#Check Bowtie:
        def test_bowtie(self):
                #print "\nSearch for bowtie"
                Search = search("bowtie")
                #print Search
                self.assertIsNotNone(Search)

        #Check Java
        def test_java(self):
                #print "\nSearch for java"
                Search = search("java")
                #print Search
                self.assertIsNotNone(Search)

        #Check Mafft:
        def test_mafft(self):
                #print "\nSearch for mafft"
                Search = search("mafft")
                #print Search
		#if Search:
		#	print "        !! Check mafft version >= 7"
                self.assertIsNotNone(Search)
	


	#Check pandas:
	def test_pandas(self):
                #print "\ntest pandas"
                Pandas_available = False
                try:
                        import pandas
                        Pandas_available = True
                except ImportError:
                        Pandas_available = False
		
		self.assertTrue(Pandas_available)
		
	#Check matplotlib:
	def test_matplotlib(self):
		#print "\ntest matplotlib"
		Matplotlib_available = False
		try:
			import matplotlib
			Matplotlib_available = True
		except ImportError:
			Matplotlib_available = False
		
		self.assertTrue(Matplotlib_available)

if __name__ == '__main__':
	#unittest.main()
	suite = unittest.TestLoader().loadTestsFromTestCase(TestConfigApytram)
	unittest.TextTestRunner(verbosity=2).run(suite)
