#!/usr/bin/python
# coding: utf-8

# File: Trinity.py
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

import logging
import subprocess
import re

def get_version():
    command = ["Trinity", "--version"]
    p = subprocess.Popen(command,
                        stdout=subprocess.PIPE,
                        stderr=subprocess.PIPE)
    (out, err) = p.communicate()
    if not err:
        version = out.split("\n")[0]
        try:
            version = map(int,version.split("v")[-1].split("."))
        except:
            version = [0]
    else:
        version = [0]
    return(version)

class Trinity(object):
    """Define an object to launch Trinity"""
    def __init__(self, OutputFile, single="", left="", right="", longreads=""):
        self.logger = logging.getLogger('apytram.lib.Trinity')
        self.seqType = "fa"
        self.max_memory = 1
        self.CPU = 1
        self.SingleInputFile = single
        self.OutputFile = OutputFile
        self.RightInputFile = right
        self.LeftInputFile = left
        self.LongReads = longreads
        self.MinLength = 200
        self.RunAsPaired = False
        self.FullCleanup = False
        self.NoPathMerging = False
        self.NoNormalizeReads = False
        self.SS_lib_type = ""

    def launch(self):
        command = ["Trinity", "--no_version_check",
                   "--seqType", self.seqType,
                   "--output", self.OutputFile,
                   "--CPU", str(self.CPU),
                   "--max_memory", str(int(self.max_memory))+"G"]

        if self.FullCleanup:
            command.append("--full_cleanup")

        if self.MinLength != 200:
            command.extend(["--min_contig_length", str(self.MinLength)])

        if self.SingleInputFile:
            command.extend(["--single", self.SingleInputFile])

        if self.RightInputFile:
            command.extend(["--right", self.RightInputFile])

        if self.LeftInputFile:
            command.extend(["--left", self.LeftInputFile])

        if self.LongReads:
            command.extend(["--long_reads", self.LongReads])

        if self.RunAsPaired:
            command.append("--run_as_paired")

        if self.NoPathMerging:
            command.append("--no_path_merging")

        if self.NoNormalizeReads:
            Trinity_version = get_version()
            if Trinity_version[0] > 2 :
                command.append("--no_normalize_reads")
            elif Trinity_version[0] == 2 and Trinity_version[1] >= 3:
                command.append("--no_normalize_reads")
            else:
                self.logger.warning("This version of Trinity has not the option --no_normalize_reads")

        if self.SS_lib_type in ["FR", "RF", "F", "R"]:
            command.extend(["--SS_lib_type", self.SS_lib_type])

        self.logger.debug(" ".join(command))
        p = subprocess.Popen(command,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        (out, err) = p.communicate()
        returncode = p.poll()

        if err:
            self.logger.error(
                 "Unexpected error when we launch Trinity:\n")
            self.logger.error(
                 " ".join(command))
            self.logger.debug(
                 "[...]\n"+"\n".join(out.strip().split("\n")[-20:]))
            self.logger.error(err)

        return (out, err, returncode)
