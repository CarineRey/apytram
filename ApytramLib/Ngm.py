#!/usr/bin/python
# coding: utf-8

# File: Ngm.py
# Created by: Carine Rey
# Created on: Feb 2017
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

import os
import logging
import subprocess


class Ngm(object):
    """Define an object to create a ngm instance"""
    def __init__(self, reference, query, sam="", output_readnames="", output_fasta = ""):
        self.logger = logging.getLogger('apytram.lib.Ngm.ngm')
        self.reference = reference
        self.query = query


        self.sensitivity = None      # --threads
        self.threads = 1          # --sensitivity
        self.min_identity = None     #-i/--min-identity
        self.min_residues = None     #-R/--min-residues
        self.min_mq = None           #-Q/--min-mq
        self.fast_pairing = False # --fast-pairing

        self.kmer      = None     #--kmer
        self.kmer_skip = None     #--kmer-skip
        self.kmer_min  = None     #--kmer-min



        self.no_unal = True        #--no-unal
        self.sam = sam             #--output
        self.output_fasta = output_fasta
        self.output_readnames = output_readnames


    def launch(self):
        command = ["ngm", "-r", self.reference, "-q", self.query, "--no-progress", "--skip-save"]

        if self.threads != 1:
            command.extend(["-t", str(self.threads)])
        if self.sensitivity != None:
            command.extend(["--sensitivity", str(self.sensitivity)])
        if self.min_identity != None:
            command.extend(["--min-identity", str(self.min_identity)])
        if self.min_residues != None:
            command.extend(["--min-residues", str(self.min_residues)])
        if self.min_mq != None:
            command.extend(["--min-mq", str(self.min_mq)])
        if self.kmer != None:
            command.extend(["--kmer", str(self.kmer)])
        if self.kmer_skip != None:
            command.extend(["--kmer-skip", str(self.kmer_skip)])
        if self.kmer_min != None:
            command.extend(["--kmer-min", str(self.kmer_min)])
        if self.no_unal:
            command.extend(["--no-unal"])


        if self.sam:
            command.extend(["--output", self.sam])
            self.logger.debug(" ".join(command))
            p = subprocess.Popen(command,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
            (out, err) = p.communicate()

        elif self.output_fasta and self.output_readnames:
            self.logger.debug(" ".join(command) + """ | awk '{OFS="\\t"; gsub("XR:i:","", $18) ;   if (( $1 !~ /^@/ ) && (($18+0) > 30)) {print ">" $1 "\\n" $10 > "%s"; print $1} }""" %self.output_fasta)
            # mapping via ngm
            p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # sam to fasta
            p1 = subprocess.Popen(["awk", """{OFS="\\t"; gsub("XR:i:","", $18) ;  if (( $1 !~ /^@/ ) && (($18+0) > 30))  {print ">" $1 "\\n" $10 > "%s"; print $1} }""" %(self.output_fasta) ], stdin=p.stdout, stdout = subprocess.PIPE)
            with open(self.output_readnames, 'w') as OUTPUTFILE:
                p2 = subprocess.Popen(["sort", "-u"], stdin=p1.stdout, stdout = OUTPUTFILE)

            p1.wait()
            p1.stdout.close()
            p2.wait()
            (out, err) = p.communicate()

        elif self.output_readnames:
            self.logger.debug(" ".join(command) + """ | cut -f 1 | grep -v -e "^@" """)
            p = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            p1 = subprocess.Popen(["cut", "-f", "1"], stdin=p.stdout, stdout = subprocess.PIPE)
            with open(self.output_readnames, 'w') as OUTPUTFILE:
                p2 = subprocess.Popen(["grep", "-v", "^@"], stdin=p1.stdout, stdout = OUTPUTFILE)

            p1.wait()
            p1.stdout.close()
            p2.wait()
            (out, err) = p.communicate()

        if os.path.isfile(self.output_readnames) and (os.stat(self.output_readnames).st_size == 0):
            if err:
                self.logger.error(err)
        return (out, err)

