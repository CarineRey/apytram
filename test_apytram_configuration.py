#!/usr/bin/python
# coding: utf-8

# File: test_configuration.py
# Created by: Carine Rey
# Created on: Nov 2015
#
#
# Copyright or © or Copr. Carine Rey
# This software is a computer program whose purpose is to assembly
# sequences from RNA-Seq data (paired-end or single-end) using one or
# more reference homologous sequences.
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

import unittest
import os
import re
import sys

import ApytramLib

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
    def test_python2_version(self):
        ver = sys.version[:3]
        self.assertTrue(ver[0]=="2")
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

    #Check Transdecoder
    def test_Transdecoder(self):
            #print "\nSearch for Transdecoder"
            Search1 = search("TransDecoder.LongOrfs")
            Search2 = search("TransDecoder.Predict")
            #print Search
            self.assertIsNotNone(Search1)
            self.assertIsNotNone(Search2)

    def test_version_Trinity(self):
            if search("Trinity"):
                Trinity_version = ApytramLib.Trinity.get_version()
                if Trinity_version[0] >= 2:
                    v = True
                    if not (Trinity_version[0] == 2 & Trinity_version[1] >= 3):
                        print "\tTrinity version must be >= v2.3 (current v%s)" %(".".join(map(str,Trinity_version)))
                else:
                    v = False
                    
                self.assertTrue(v)

    #Check samtools:
    def test_samtools(self):
        Search = search("samtools")
        if Search:
            pass
        elif search("Trinity"):
            Trinity_version = ApytramLib.Trinity.get_version()
            if Trinity_version[0] < 2:
                pass
            elif (Trinity_version[0] == 2 & Trinity_version[1] <= 2):
                pass
            else:
                print("\tSamtools not needed")
                Search = "ok"
        self.assertIsNotNone(Search)

    #Check Bowtie2:
    def test_bowtie2(self):
            #print "\nSearch for bowtie2"
            Search = search("bowtie2")
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
        #   print "        !! Check mafft version >= 7"
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
    ret = not unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()
    sys.exit(ret)
