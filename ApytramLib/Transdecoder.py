#!/usr/bin/python
# coding: utf-8

# File: Transdecoder.py
# Created by: Carine Rey
# Created on: Nov 2017
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

import logging
import subprocess
import re

def get_version():
    command = ["TransDecoder.LongOrfs", "--version"]
    p = subprocess.Popen(command,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
    (out, err) = p.communicate()
    if not err:
        version = out.split(" ")[1]
        try:
            version = map(int,version.split("."))
        except:
            version = [0]
    else:
        version = [0]
    return(version)

class TransDecoder(object):
    """Define an object to launch TransDecoder"""
    def __init__(self, FastaFile, min_prot_length=None, cpu=1, single_best_orf=False, TmpDir=None):
        self.logger = logging.getLogger('apytram.lib.TransDecoder')
        self.fasta = FastaFile
        self.min_prot_length = min_prot_length
        self.retain_long_orfs = min_prot_length
        self.single_best_orf = single_best_orf
        self.cpu = cpu
        self.TmpDir = TmpDir

    def launch(self):
        command1 = ["TransDecoder.LongOrfs", "-t", self.fasta]
        command2 = ["TransDecoder.Predict", "-t", self.fasta]

        command2.extend(["--no_refine_starts"])
        if self.min_prot_length:
            command1.extend(["-m", str(int(self.min_prot_length/3 - 1))])
        if self.min_prot_length:
            command2.extend(["--retain_long_orfs_mode", str("strict")])
            command2.extend(["--retain_long_orfs_length", str(int(self.min_prot_length))])
        if self.single_best_orf:
            command2.append("--single_best_orf")

        self.logger.debug(" ".join(command1))
        p1 = subprocess.Popen(command1,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE,
                             cwd=self.TmpDir)
        (out1, err1) = p1.communicate()
        returncode1 = p1.poll()

        returncode2 = 1
        if returncode1:
            self.logger.error(
                 "Unexpected error when we launch TransDecoder.LongOrfs:\n")
            self.logger.error(
                 " ".join(command1))
            self.logger.debug(
                 "[...]\n"+"\n".join(out1.strip().split("\n")[-20:]))
            self.logger.error(err1)

        else:
            p2 = subprocess.Popen(command2,
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                cwd=self.TmpDir)
            (out2, err2) = p2.communicate()
            returncode2 = p2.poll()

            if returncode2:
                self.logger.error(
                    "Unexpected error when we launch TransDecoder.Predicts:\n")
                self.logger.error(
                    " ".join(command2))
                self.logger.debug(
                    "[...]\n"+"\n".join(out2.strip().split("\n")[-20:]))
                self.logger.error(err2)

        return (returncode1+returncode2)
