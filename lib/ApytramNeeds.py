import os
import sys
import subprocess

def fastq2fasta(FastqFile,FastaFile):
    ExitCode = 1
    command = """awk 'NR%%4==1||NR%%4==2' < %s | tr "@" ">" > %s """ %(FastqFile, FastaFile)
    os.system(command)
    return ExitCode

def add_paired_read_names(File):
    "Add paired read name to a read name list"
    if os.path.isfile(File):
        command1 = """awk '{ print $0; if (match($0,"1$")) sub("1$",2,$0); else if (match($0,"2$")) sub("2$",1,$0); print $0}' %s  | sort -u > %s """ %(File, File+".paired")
        command2 = "mv %s %s" %(File+".paired", File)
        os.system(command1)
        os.system(command2)
    return 0
