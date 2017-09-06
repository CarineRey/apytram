#!/usr/bin/python
# coding: utf-8

# File: ApytramNeeds.py
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


import os
import re
import shutil
import string
import numpy as np
import sys
import subprocess
import pandas
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#matplotlib.style.use('ggplot')



# Class

#### Sequence class

class Sequence(object):
    def __init__(self,):
        self.Name = ""
        self.OldName = ""

        self.Sequence = ""

        ### homology score with the best reference
        #t = target
        #q = query
        #i = id
        #l = length
        #a = aligned part
        #b = begin
        #e = end
        self.ti =    ""
        self.ql =    ""
        self.tal =   ""
        self.tl =    ""
        self.tab =   ""
        self.tae =   ""
        self.score = ""
        self.pi =    ""
        self.qab =   ""
        self.qae =   ""

        ### Caracteristics
        self.BestSequence = ""
        self.Complement = False

    def __str__(self):
        return(">" + self.Name + "\n" + '\n'.join(self.Sequence[i:i+60] for i in range(0, len(self.Sequence), 60)) + "\n")

    def __len__(self):
        return len(self.Sequence)

    def copy(self):
        ns = Sequence()

        for (att,value) in self.__dict__.items():
            ns.__dict__[att] = value

        return(ns)

    def add_attribute_from_exonerate(self, ExonerateListLine):
        #TrinityExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%tal\t%tl\t%tab\t%tae\t%s\t%pi\t%qab\t%qae\n"
        #"%ti\t%qi\t%ql\t%tal\t%tl\t%tab\t%tae\t%s\t%pi\t%qab\t%qae\n"
        self.ti     = ExonerateListLine[0]
        self.Name   = ExonerateListLine[1]
        self.ql     = float(ExonerateListLine[2])
        self.tal    = float(ExonerateListLine[3])
        self.tl     = float(ExonerateListLine[4])
        self.tab    = float(ExonerateListLine[5])
        self.tae    = float(ExonerateListLine[6])
        self.score  = float(ExonerateListLine[7])
        self.pi     = float(ExonerateListLine[8])
        self.qab    = float(ExonerateListLine[9])
        self.qae    = float(ExonerateListLine[10])

    def complete(self,newinfo_sequence):
        for (att,value) in self.__dict__.items():
            if not value and newinfo_sequence.__dict__[att]:
                self.__dict__[att] = newinfo_sequence.__dict__[att]



#### Fasta class

# import os,re
# f = Fasta()
# f.read_fasta(FastaFilename = "example/multi_ref.fasta")
# fi = f.filter_fasta("ENSMUSP00000002708_G01455")
# print fi
# d = {"ENSMUSP00000002708_G01455" : "totot"}
# fd = fi.rename_fasta(d)
# print fd
# fd.write_fasta("test")

class Fasta(object):
    def __init__(self):
        self.Sequences = []
        self.Names = []
        self.ReferenceSequences = []

    def __str__(self):
        string = []
        for s in self.Sequences:
            string.extend([str(s)])
        return("".join(string))

    def append(self,new_sequence):
        assert isinstance(new_sequence, Sequence), "Sequence must belong to the Sequence class"
        if new_sequence.Complement and new_sequence.Sequence:
            new_sequence.Sequence = complement(new_sequence.Sequence)
            new_sequence.Complement = False

        self.Sequences.append(new_sequence)
        self.Names.append(new_sequence.Name)

    def read_fasta(self, FastaFilename = "" , String = ""):
        if String:
            Fasta = String.strip().split("\n")
        elif os.path.isfile(FastaFilename):
            File = open(FastaFilename,"r")
            Fasta = File.read().strip().split("\n")
            File.close()
        else:
            Fasta = []

        name = ""
        sequence = ""
        sequence_list = []

        for line in Fasta + [">"]:
            if re.match(">",line):
                # This is a new sequence write the previous sequence if it exists
                if sequence_list:
                    new_sequence = Sequence()
                    new_sequence.Name = name
                    new_sequence.Sequence = "".join(sequence_list)
                    self.append(new_sequence)
                    sequence_list = []

                name = line.split()[0][1:] # remove the >

            elif name != "":
                sequence_list.append(line)
            else:
                pass

    def filter_fasta(self, SelectedNames):
        FilteredFasta = Fasta()
        for s in self.Sequences:
            if s.Name in SelectedNames:
                FilteredFasta.append(s)
        return FilteredFasta

    def dealign_fasta(self):
        DealignedFasta = Fasta()
        for s in self.Sequences:
            s.Sequence = s.Sequence.replace("-", "")
            DealignedFasta.append(s)
        return DealignedFasta

    def isalign(self):
        len_s = []
        gap_s = []
        isalign = False
        for s in self.Sequences:
            len_s.append(len(s))
            if len_s and len(s) != len_s[0]:
                break
            elif re.search("-",s.Sequence):
                gap_s.append(True)
        if len(set(len_s)) == 1 and gap_s:
            isalign = True
        return isalign

    def complete_fasta(self,NewInfoFasta):
        assert isinstance(NewInfoFasta, Fasta), "NewInfoFasta must belong to the Fasta class"
        for i in range(len(self.Sequences)):
            if self.Sequences[i].Name in NewInfoFasta.Names:
                index = NewInfoFasta.Names.index(self.Sequences[i].Name)
                self.Sequences[i].complete(NewInfoFasta.Sequences[index])


    def rename_fasta(self, OldNamesNewNamesDict):
        assert isinstance(OldNamesNewNamesDict, dict), "OldNamesNewNamesDict must be a dict"
        RenamedFasta = Fasta()
        for s in self.Sequences:
            if s.Name in OldNamesNewNamesDict.keys():
                ns = s.copy()
                ns.Name = OldNamesNewNamesDict[s.Name]
                ns.OldName = s.Name
                RenamedFasta.append(ns)
        return RenamedFasta


    def write_fasta(self, OutFastaFile):
        # Write all sequences in the file
        write_in_file(str(self), OutFastaFile)


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

def end(exit_code, TmpDirName, keep_tmp = False, logger=""):
    "Functions to end apytram in removing temporary directory"
    ### Remove tempdir if the option --tmp have not been use
    if logger:
        logger.debug("Remove the temporary directory")
    #Remove the temporary directory :
    if not keep_tmp and "apytram" in TmpDirName:
        shutil.rmtree(TmpDirName)
    sys.exit(exit_code)


def set_directory_from_prefix(Prefix,DirType="",logger=""):
        ### Set up the output directory
        if Prefix:
            DirName = os.path.dirname(Prefix)
            if os.path.isdir(DirName):
                if DirType and logger:
                    logger.info("The %s directory %s exists" %(DirType,DirName))
            elif DirName: # if DirName is not a empty string we create the directory
                if DirType and logger:
                    logger.info("The %s directory %s does not exist, it will be created" %(DirType,DirName))
                os.makedirs(DirName)


def fastq2fasta(FastqFiles,FastaFile):
    ExitCode = 1
    command1 = """cat %s """ %(FastqFiles)
    command2 = ["awk", """NR%4==1||NR%4==2"""]
    command3 = ["tr", "@", ">"]
    with open(FastaFile, 'w') as OutFile:
        p1 = subprocess.Popen(command1.split(), stdout=subprocess.PIPE)
        p2 = subprocess.Popen(command2, stdin=p1.stdout, stdout=subprocess.PIPE)
        p3 = subprocess.Popen(command3, stdin=p2.stdout, stdout=OutFile)
        p1.stdout.close()
        p2.stdout.close()
        out, err = p3.communicate()
        p1.wait()
        p2.wait()
    return (out, err)

def cat_fasta(FastaFiles,CatFastaFile):
    command = "cat %s" %(FastaFiles)
    with open(CatFastaFile, 'w') as OutFile:
        p = subprocess.Popen(command.split(), stdout = OutFile)
        (out, err) = p.communicate()
    return (out, err)

def complement(Sequence_str):
    intab = "ABCDGHMNRSTUVWXYabcdghmnrstuvwxy"
    outtab = "TVGHCDKNYSAABWXRtvghcdknysaabwxr"
    trantab = string.maketrans(intab, outtab)
    # Reverse
    Reverse = Sequence_str.replace("\n","")[::-1]
    # Complement
    Complement = Reverse.translate(trantab)
    return Complement

def write_in_file(String,Filename,mode = "w"):
    if mode in ["w","a"]:
        with open(Filename,mode) as File:
            File.write(String)




def read_clstr(File):
    res = {}
    if os.path.isfile(File):
        File = open(File,"r")
        Fasta = File.read().strip().split("\n")
        File.close()
    else:
        return res

    name = ""
    sequence = ""
    sequence_list = []

    for line in Fasta + [">Cluster"]:
        if re.match(">Cluster",line):
            # This is a new clustr write the previous sequence if it exists
            if sequence_list:
                res[name] = sequence_list
                sequence_list = []

            name = line[1:].replace("Cluster ","Cluster_") # remove the >

        elif name != "":
            seq = line.split(">")[1].split("...")[0]
            rep = line.split(" ")[-1] == "*"
            if not rep:
                rev = line.split(" ")[-1].split("/")[0] == "-"
            else:
                rev = False
            sequence_list.append((seq, rep, rev))
        else:
            pass
    return res



def add_paired_read_names(File, NewFile, logger = ""):
    "Add paired read names to a read name list"
    if os.path.isfile(File):
        command1 = ["awk", """{ print $0; if (match($0,"1$")) sub("1$",2,$0); else if (match($0,"2$")) sub("2$",1,$0); print $0}""", File]
        command2 = "sort -u -o %s" %(NewFile)
        p1 = subprocess.Popen(command1, stdout = subprocess.PIPE)
        p2 = subprocess.Popen(command2.split(), stdin = p1.stdout,
                              stdout = subprocess.PIPE )
        p1.stdout.close()
        (out, err) = p2.communicate()
        p1.wait()

        if err:
            logger.error(err)
    elif logger:
        logger.error("%s is not a file" %(File))

def remove_duplicated_read_names(File, NewFile, logger = ""):
    "remove duplicated read names"
    if os.path.isfile(File):
        command = """sort -u -o %s %s""" %(File, NewFile)
        p = subprocess.Popen(command.split(), stdout = subprocess.PIPE)
        (out, err) = p.communicate()
        if err:
            logger.error(err)
    elif logger:
        logger.error("File %s is not a valid path" %File)


def get_free_space(DirPath):
    if not DirPath:
        DirPath = "."
    free = 0
    if os.path.isdir(DirPath):
        st = os.statvfs(DirPath)
        free = (st.f_bavail * st.f_frsize)
    return free

def get_size(Files):
    size = 0
    if not type(Files) == type([]):
        Files = [Files]
    for f in Files:
        if os.path.isfile(f):
            size += os.stat(f).st_size
    return size

def count_lines(Files):
    Number = 0
    if type(Files) != type([]):
        Files = [Files]
    for File in Files:
        if os.path.isfile(File):
            with open(File, 'r') as InFile:
                #Exit = subprocess.check_output(["wc", "-l", File])
                command = "wc -l"
                p = subprocess.Popen(command.split(),
                                 stdin=InFile,
                                 stdout=subprocess.PIPE)
                (out, err) = p.communicate()
                Number += int(out.strip())
    return Number

def count_sequences(File):
    out = "0\n"
    if os.path.isfile(File):
        #Exit = subprocess.check_output(["grep \"^>\" %s | wc -l" %File], shell=True)
        #Number = Exit.replace("\n","")
        command = """grep -c ^> %s""" %(File)
        p = subprocess.Popen(command.split(),
                             stdout=subprocess.PIPE)
        (out, err) = p.communicate()
    return int(out.strip())

def split_readnames_in_right_left(File,Right,Left):
    if os.path.isfile(File):
        command1 = ["grep", "1$", File]
        command2 = ["grep", "2$", File]
        with open(Right, 'w') as OutFileR:
            p1 = subprocess.Popen(command1, stdout = OutFileR)
            (_, err1) = p1.communicate()
        with open(Left, 'w') as OutFileL:
            p2 = subprocess.Popen(command2, stdout = OutFileL)
            (_, err2) = p2.communicate()

        return ("","")
    else:
        return("","File %s is not a valid path" %File)

def check_paired_data(FastaFile):
    BadReadName = ""
    if os.path.isfile(FastaFile):
        #Check First sequence name end with 1/ or 2/
        command1=["grep", "^>",FastaFile, "-m", "100"]
        command21=["tail", "-n" , "300",  FastaFile]
        command22=["grep",  "-m",  "100",  """^>"""]

        p1 = subprocess.Popen(command1, stdout = subprocess.PIPE)
        (out1, err1) = p1.communicate()

        with open(os.devnull, 'w')  as FNULL:
            p21 = subprocess.Popen(command21, stdout = subprocess.PIPE,
                            stderr = FNULL)
        p22 = subprocess.Popen(command22, stdin=p21.stdout,
                               stdout = subprocess.PIPE)

        p21.stdout.close()
        (out2, err2) = p22.communicate()
        p21.wait()

        #Exit = subprocess.check_output(["grep", "^>",FastaFile, "-m", "100"])
        #Exit2 = subprocess.check_output(["tac %s | grep -m 100 \"^>\"" %FastaFile], shell=True)
        ReadNames = (out1+out2).strip().split("\n")
        while (not BadReadName) and ReadNames:
            ReadName = ReadNames.pop()
            if not re.search("[12]$",ReadName):
                BadReadName = ReadName[1:]
    return BadReadName

def are_identical(File1, File2):
    Identical = False
    if os.path.isfile(File1) and  os.path.isfile(File2):
        with open(os.devnull, 'w')  as FNULL:
            diff = subprocess.call(["diff",File1,File2,"--brief"],
                                       stdout=FNULL,
                                       stderr=FNULL)
        if not diff:
            Identical = True
    return Identical

def number_new_reads(FileOld, FileNew, nb_intial=0):
    Nb = 0
    if os.path.isfile(FileNew) and  os.path.isfile(FileOld):
        command1 = "join -v 2 %s %s" %(FileOld, FileNew)
        command2 = "wc -l"
        p1 = subprocess.Popen(command1.split(),
                             stdout=subprocess.PIPE)
        p2 = subprocess.Popen(command2.split(),
                             stdin=p1.stdout,
                             stdout=subprocess.PIPE)
        # Allow p2 to receive a SIGPIPE if p1 exits.
        p1.stdout.close()
        (out, err) = p2.communicate()
        p1.wait()

        if out:
           Nb = int(out)
        else:
            print err
    elif nb_intial:
        Nb = nb_intial
    return Nb

def new_reads(FileOld, FileNew, DiffFile,f=""):
    if os.path.isfile(FileNew) and os.path.isfile(FileOld):
        command1 = "join -v 2 %s %s" %(FileOld, FileNew)
        with open(DiffFile, "w") as OUTPUTFILE:
            p1 = subprocess.Popen(command1.split(),
                             stdout=OUTPUTFILE)
        (out, err) = p1.communicate()
        return(count_lines(DiffFile), DiffFile)
    elif os.path.isfile(FileNew):
        return(count_lines(FileNew), FileNew)
    else:
        return(0, "")

def common_reads(File1, File2, CommonFile):
    if os.path.isfile(File1) and os.path.isfile(File2):
        command1 = "join %s %s" %(File1, File2)
        with open(CommonFile, "w") as OUTPUTFILE:
            p1 = subprocess.Popen(command1.split(),
                             stdout=OUTPUTFILE)
        (out, err) = p1.communicate()
        return(count_lines(CommonFile))
    else:
        return(0)


def get_seq(id, DB):
    try:
        return DB.get_raw(id)
    except KeyError:
        return ""

def retrieve_reads_from_index(DB, ReadsNamesFile, OutputFile):
    reads_file = open(ReadsNamesFile, "r")
    reads = reads_file.read().strip().split("\n")
    reads_file.close()
    res=[]
    for r in reads:
        res.append(get_seq(r, DB))
    with open(OutputFile, "w") as FILE:
        FILE.write("".join(res))
    return(len(res))

def get_read_from_cluster(rep, DB):
    try:
        cluster = DB.get_raw(rep)
        return ([rep] + cluster.split("\n")[1].split(";"))
    except KeyError:
        return [rep]

def retrieve_reads_from_cluster(DB, ReadsNamesFile, OutputFile):
    reads_file = open(ReadsNamesFile, "r")
    reads = reads_file.read().strip().split("\n")
    reads_file.close()
    res=[]
    for r in reads:
        res.extend(get_read_from_cluster(r, DB))
    with open(OutputFile, "w") as FILE:
        FILE.write("\n".join(res))
    return(len(res))

def check_almost_identical_exonerate_results(ExonerateResult):
    "Return True if all hit have a hit with 99% > id and a len = 98% len query "
    Result = False
    i = 0
    List = ExonerateResult.strip().split("\n")
    for line in List:
        ListLine = line.split("\t")
        if len(ListLine) > 7:
            #ExonerateProcess.Ryo = "%ti\t%qi\t%ql\t%qal\t%tal\t%tl\t%pi\t\n"
            pi = float(ListLine[6])
            ql = float(ListLine[2])
            tl = float(ListLine[5])
            pl = min( ql, tl) / max( ql, tl)
            if (pl >= 0.99) and (pi > 98) :
                i += 1

    if len(List) == i:
        Result = True
    return Result

def write_stats(df_list, Output):
    df = pandas.concat(df_list)
    df.to_csv(Output, na_rep = "0.0",index=False)

def create_plot(TimeStatsDictList, IterStatsDictList, SpeciesListNames, OutPrefixName):
    with PdfPages("%s.stats.pdf" %(OutPrefixName)) as pdf:
        for i in range(len(TimeStatsDictList)):

            TimeStatsDict = TimeStatsDictList[i]

            # remove keys without value
            for iter_i in TimeStatsDict.keys():
                for att in TimeStatsDict[iter_i].keys():
                    if TimeStatsDict[iter_i][att] == 0:
                        del TimeStatsDict[iter_i][att]
            IterStatsDict = IterStatsDictList[i]
            species = SpeciesListNames[i]

            df_time = pandas.DataFrame(TimeStatsDict).T
            df_iter = pandas.DataFrame(IterStatsDict).T


            ### Plot 1 ###
            df_iter.plot(subplots=True,
                        grid = True,
                        legend = "best",
                        figsize=(10, 20),
                        style = "-o",
                        title = 'Evolution of some values - %s' %(species))
            plt.xlabel('Iteration')
            pdf.savefig()
            plt.close()

            ### Plot 2 ###
            df_time.plot(kind ='bar',
                        grid = True,
                        figsize=(10, 10),
                        stacked=True,
                        legend = "best")
            plt.ylabel('seconds')
            plt.xlabel('Iteration')
            plt.title('Distribution of the execution time - %s' %(species))
            pdf.savefig()
            plt.close()


def calculate_coverage(AlignmentString, RefSeq):
    assert AlignmentString, "AlignmentString must not be empty"
    assert isinstance(AlignmentString,str),  "AlignmentString must be str"
    assert isinstance(RefSeq,list),  "RefSeq must be list"
    # The first sequence must be the reference
    #Alignment reading
    Alignment = Fasta()
    Alignment.read_fasta(String = AlignmentString)

    #Build a reference consensus sequences

    Ref_Ali = Alignment.filter_fasta(RefSeq)
    ConsensusList = []
    AliLength = len(Ref_Ali.Sequences[0].Sequence)
    for i in range(AliLength):
        # get all nucl for reference for each position
        pos = {seq.Sequence[i] for seq in Ref_Ali.Sequences}
        if "-" in pos:
            pos.remove("-")
        pos = list(pos)
        if len(pos) == 1:
            ConsensusList.append(pos[0])
        elif len(pos) == 0:
            ConsensusList.append("-")
        else:
            ConsensusList.append("0")

    # Attribute a code for each position
    ## 0 Gap
    ## 1 All references identical
    ## 2 ambigous in reference
    ## 3 contig identical references
    ## 4 contig ok but ambigous reference
    ## 5 present in contig but not in reference
    ## 6 different in contig than in reference

    DicPlotCov = {}
    StrictRefCov = [0]*AliLength
    LargeRefCov = [0]*AliLength

    for Sequence in Alignment.Sequences:
        Ref = False
        Name = Sequence.Name
        seq = Sequence.Sequence
        DicPlotCov[Name] = {}

        if Name in RefSeq:
            Ref = True

        for i in range(AliLength):

            if Ref:
                if ConsensusList[i] == "0":
                    DicPlotCov[Name][i] = 2
                elif seq[i] == "-":
                    DicPlotCov[Name][i] = 0
                else:
                    DicPlotCov[Name][i] = 1
            else:
                if seq[i] == "-":
                    DicPlotCov[Name][i] = 0
                elif ConsensusList[i] == "-":
                    DicPlotCov[Name][i] = 5
                    LargeRefCov[i] = 1
                elif ConsensusList[i] == seq[i]:
                    DicPlotCov[Name][i] = 3
                    StrictRefCov[i] = 1
                    LargeRefCov[i] = 1
                elif ConsensusList[i] == "0":
                    DicPlotCov[Name][i] = 4
                    StrictRefCov[i] = 1
                    LargeRefCov[i] = 1
                else:
                    DicPlotCov[Name][i] = 6
                    StrictRefCov[i] = 1
                    LargeRefCov[i] = 1

    # Get length of the consensus
    ConsensusLength = len([ p for p in ConsensusList if p != "-" ])

    # We count the number of base position in the alignment
    SumStrictRefCov = sum(StrictRefCov)
    SumLargeRefCov = sum(LargeRefCov)

    #Percentage of positon of the reference which is reprensented in contigs
    StrictCov=float(SumStrictRefCov)/float(ConsensusLength) *100
    #Percentage of positon in contigs on the number positon of the reference
    #It can be superior to 100 if the contigs are longer than the reference
    LargeCov=float(SumLargeRefCov)/float(ConsensusLength) *100
    return (StrictCov, LargeCov, DicPlotCov)

def create_plot_ali(DictPlotCov, OutPrefixName):
    ### Plot Ali ###
    "Create a png file containing a representation of an alignement. The first seqeunce must be the reference."
    df = pandas.DataFrame(DictPlotCov).T
    fig, ax = plt.subplots()
    ax.grid(False)
    cmap, norm = matplotlib.colors.from_levels_and_colors([0,0.5, 1.5, 2.5,3.5,4.5,5.5,6.5],
                                                          ["White","Blue","DarkViolet",'Darkgreen',"DarkViolet","Orange",'Darkred'])
    heatmap = ax.pcolor(df ,
                        edgecolors="white",  # put black lines between squares in heatmap
                        cmap=cmap,
                        norm=norm)
    # Format
    fig = plt.gcf()
    Width = 3 + float(df.shape[1])/12
    Heigth =  1+float(df.shape[0])/2

    # Maximum figure size (pixels)
    if Width > 32000:
        Width = 32000
    if Heigth > 32000:
        Heigth = 32000

    fig.set_size_inches(Width, Heigth)
    # put the major ticks at the middle of each cell
    plt.yticks(np.arange(len(df.index)) + 0.5, df.index, size = 12)
    ax.xaxis.tick_top()
    plt.xticks(np.arange(len(df.columns)) + 0.5, df.columns, rotation=90, size= 6)
    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.tick_params(bottom='off', top='off', left='off', right='off')
    plt.xlabel('Position')
    fig.tight_layout()
    fig.savefig("%s.ali.png" % OutPrefixName)
    plt.close()
