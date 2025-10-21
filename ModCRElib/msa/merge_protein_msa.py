import os, sys, re
from collections import Counter
import configparser
import itertools
import numpy
import optparse
import subprocess
import time
import random


# Get scripts path (i.e. ".") #
exe_path = os.path.abspath(os.path.dirname(__file__))
if os.path.exists(os.path.join(exe_path,"..","ModCRElib")):
   scripts_path = os.path.join(exe_path,"..")
elif os.path.exists(os.path.join(exe_path,"..","..","ModCRElib")):
   scripts_path = os.path.join(exe_path,"..","..")
elif os.path.exists(os.path.join(exe_path,"..","..","..","ModCRElib")):
   scripts_path = os.path.join(exe_path,"..","..","..")
else:
   scripts_path = os.path.join(exe_path)

config_path  = os.path.join(scripts_path,"ModCRElib","configure")


# Append scripts path to python path #
sys.path.append(scripts_path)

# Read configuration file #
config = configparser.ConfigParser()
config_file = os.path.join(config_path, "config.ini")
config.read(config_file)

# Imports my functions #
from ModCRElib.beans import functions

# Imports jbonet's module #
from SBILib.data import aminoacids1to3, aminoacids_polarity_boolean, nitrogenous_bases,dna_complementary
from SBILib.structure import PDB

# Import my modules #
from ModCRElib.msa import pwm_pbm as PWM

# Defined maximum size of nMSA, minimum must be 1000
try:
  maxsize= int(config.get("Parameters", "max_sequences_in_MSA"))
except:
  maxsize=5000
if maxsize < 1000: maxsize=1000


#-------------#
# Functions   #
#-------------#

def parse_options():
    """
    This function parses the command line arguments and returns an optparse
    object.
    
    """

    parser = optparse.OptionParser("python merge_protein_msa.py -i input_main_alignment -j input_list_msa -w weight [-o output_file --dummy dummy_dir ]  ")

    parser.add_option("--dummy", default="/tmp/", action="store", type="string", dest="dummy_dir", help="Dummy directory (default = /tmp/)", metavar="{directory}")
    parser.add_option("-i", "--alignment", action="store", default=None, type="string", dest="input_main_alignment", help="Alignment of modelled TF sequences associated with structural models. Format FastA, it can be produced by Superimposition and must not have internal gaps ", metavar="{filename}")
    parser.add_option("-j", "--MSA_list", action="store", default=None, type="string", dest="input_list_msa", help="List of MSAs (format: code  path/file  weigth). Code is used to identify the corresponding sequence in the alignment file", metavar="{filename}")
    parser.add_option("-m", "--maxsize", action="store", default=None, type="int", dest="maxsize", help="Maximum number of sequences in a MSA file (default uses max_sequences_in_MSA of configuration) ", metavar="{filename}")
    parser.add_option("-o", action="store", default="output_msa", type="string", dest="output_file", help="Output MSA (default = 'output_msa')", metavar="{filename}")
    parser.add_option("-v", "--verbose", default=False, action="store_true", dest="verbose", help="Verbose mode (default = False)", metavar="{boolean}")
  
 
    (options, args) = parser.parse_args()

    if options.input_main_alignment is None or options.input_list_msa is None :
        parser.error("missing arguments: type option \"-h\" for help")

    return options

def get_offsets(fasta_file):
    offsets={}
    for header, sequence in functions.parse_fasta_file(fasta_file):
        #print(header,sequence,len(sequence),len(sequence.lstrip("X")),len(sequence.rstrip("X")))
        offsets.setdefault(  header,(    ( len(sequence) -len(sequence.lstrip("X")) ),(  len(sequence) -len(sequence.rstrip("X"))  )  )  ) 
    return offsets

def get_msas(msa_list):
    msas={}
    for line in functions.parse_file(msa_list):
        code,msa_file,weight = line.split()
        try:
           msa=PWM.pMSA(file_name=msa_file, motif_name=code, option="msa")
        except:
           print("Warning, not found file %s"%msa_file)
           try:
              msa=PWM.pMSA(file_name=os.path.join(os.path.dirname(msa_list),msa_file), motif_name=code, option="msa")
           except Exception as e:
              raise(e)
        msas.setdefault(code,(msa,weight))
    return msas

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    # Arguments & Options #
    options = parse_options()
    fasta_file=options.input_main_alignment
    msa_list  =options.input_list_msa
    output    =options.output_file
    verbose   =options.verbose
    dummy_dir =options.dummy_dir

    if options.maxsize is not None:
       maxsize=int(options.maxsize)

    offsets = get_offsets(fasta_file)
    #print(offsets)
    try:
      msas    = get_msas(msa_list)
    except Exception as e:
      print(e)
      exit(0)
    wsum = sum([float(msas[code][1]) for code in msas])
    merge_msa  = PWM.pMSA()
    for code in msas:
        if code not in offsets: continue
        msa,weight=msas[code]
        start,last = offsets[code]
        n = maxsize*float(weight)/wsum 
        i=0
        while (i<n):
          for seq_scr in sorted(msa.get_sequences(),key=lambda x: x[1]):
            if i<=n:
               merge_msa.add_sequence(start*"."+seq_scr[0]+last*".",seq_scr[1])
               i+=1
    merge_msa.set_motif(os.path.basename(output))
    if verbose: sys.stdout.write("-- Write meme file %s...\n"%(output.rstrip(".msa")+".meme"))
    merge_msa.write(output.rstrip(".msa")+".meme",option="meme")
    if verbose: sys.stdout.write("-- Write PWM file %s...\n"%(output.rstrip(".msa")+".pwm"))
    merge_msa.write(output.rstrip(".msa")+".pwm",option="pwm")
    if verbose: sys.stdout.write("-- Write MSA file %s...\n"%(output.rstrip(".msa")+".msa"))
    merge_msa.write(output.rstrip(".msa")+".msa",option="msa")
    if verbose: sys.stdout.write("-- Write LOGOS %s...\n"%(output.rstrip(".msa")+".logo"))
    PWM.write_protein_logo(merge_msa,output.rstrip(".msa")+".logo",dummy_dir)
    if verbose: print("Done")    
    



    
