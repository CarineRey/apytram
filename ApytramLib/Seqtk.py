#!/usr/bin/python
# coding: utf-8

# File: Seqtk.py
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


class Seqtk(object):
    """Define an object to create a seqtk instance"""
    def __init__(self, fasta="", fastq=""):
        self.logger = logging.getLogger('apytram.lib.Seqtk.seqtk')
        self.fastq = fastq
        self.fasta = fasta

    def launch_fastq2fatsa(self, output):
        command = ["seqtk", "seq", "-A", self.fastq]
        with open(output, 'w') as OUTPUTFILE:
            self.logger.debug(" ".join(command))
            p = subprocess.Popen(command, stdout = OUTPUTFILE)
            (out, err) = p.communicate()
        return (out, err)

    def launch_fasta_subseq(self, list_seqs, output):
        command = ["seqtk", "subseq",  self.fasta, list_seqs]
        with open(output, 'w') as OUTPUTFILE:
            self.logger.debug(" ".join(command))
            p = subprocess.Popen(command, stdout = OUTPUTFILE)
            (out, err) = p.communicate()
        return (out, err)

