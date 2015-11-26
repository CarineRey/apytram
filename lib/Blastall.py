import os
import sys
import logging
import subprocess

## TO DO: check if blastall is installed and is path is blastall


class Formatdb:
    """Define an object to create a local database"""
    def __init__(self, InputFile, OutputFiles):
        self.logger = logging.getLogger('Blastplus.Formatdb')
        self.logger.info('creating an instance of Formatdb')
        self.Title = ""
        self.InputFile = InputFile
        self.OutputFiles = OutputFiles
        self.LogFile = ""
        self.Protein = "F" # "F" => nucleotides "T" => Proteins
        self.IndexedDatabase = "T" # "FT"

    def launch(self):
        self.logger.info('doing something')
        command = ["formatdb","-i",self.InputFile,"-n", self.OutputFiles,
                   "-p", self.Protein, "-o" , self.IndexedDatabase]
        self.logger.info(" ".join(command))
        try:
            result = subprocess.call(command)
        except:
            os.system("echo Unexpected error: "+" ".join(command)+"\n")
        return command


class Blastall:
    """Define a object to lauch a blast on a local database"""
    def __init__(self, Program, Database, QueryFile):
        self.Program = Program
        self.Database = Database
        self.QueryFile = QueryFile
        self.OutputFile = ""
        self.Threads = 1
        self.OutFormat = 8
        self.Evalue = 10
        #-v  Number of database sequences to show one-line descriptions for (V) [Integer]
        #default = 500
        self.v = 100000000000
        #-b  Number of database sequence to show alignments for (B) [Integer]
        #default = 250
        self.b = 100000000000

    
    def launch(self,OutputFile):
        command = ["blastall","-p",self.Program,"-d", self.Database, "-i" , os.path.abspath(self.QueryFile),
                   "-o", os.path.abspath(OutputFile), "-e", str(self.Evalue), "-m", str(self.OutFormat),
                  "-a", str(self.Threads)] #, "-v", str(self.v), "-b", str(self.b) ]
        try:
            #TO DO recuperer dans des fichier ouvert open wb balise et balise Err les sorties standards
#balise=open('./OUT.txt', 'wb')
#baliseErr=open('./ERR.txt', 'wb')
#result = subprocess.call(command, stdout=balise,stderr=baliseErr)
                subprocess.call(command)
        except:
            os.system("echo Unexpected error\n")
            print " ".join(command)
        return command
    
    def get_hit_names(self,OutputFile):
        command = ["blastall","-p",self.Program,"-d", self.Database, "-i" , os.path.abspath(self.QueryFile),
                  "-e", str(self.Evalue), "-m", str(self.OutFormat),
                  "-a", str(self.Threads)] #, "-v", str(self.v), "-b", str(self.b) ]
        try:
            #TO DO recuperer dans des fichier ouvert open wb balise et balise Err les sorties standards
#balise=open('./OUT.txt', 'wb')
#baliseErr=open('./ERR.txt', 'wb')
#result = subprocess.call(command, stdout=balise,stderr=baliseErr)
            result = subprocess.check_output(command)
        except:
            os.system("echo Unexpected error when we launch blastall with :\n")
            result = ""
            print " ".join(command)
        
        #Get hit names:
        Names = [x.split("\t")[1] for x in result.strip().split("\n")]
        #Write Names in the outputfile
        try:
            Output = open(OutputFile,"w")
            Output.write("\n".join(Names))
            Output.close()
        except:
            os.system("echo Unexpected error when we tried to write in %s\n" % OutpuFile)
        
        return Names


class Fastacmd:
    """Define a object to launch fastacmd on a local database"""
    def __init__(self, Database, SequenceNamesFile, OutputFile):
        self.InputFile = SequenceNamesFile
        self.Database = Database
        self.OutputFile = OutputFile
        self.Protein = "F"

    def launch(self):
        command = ["fastacmd","-d",self.Database,"-i", self.InputFile,
                   "-p", self.Protein, "-o" , self.OutputFile ]
        print " ".join(command)
        
        try:
            subprocess.call(command)
        except:
            os.system("echo Unexpected error when we launched fastacm:\n")
            print " ".join(command)
        return command