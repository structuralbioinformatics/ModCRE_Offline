from SBILib.structure import PDB
import sys
file=sys.argv[1]

protein = PDB(file)

for x in protein.proteins:
    begin = x.first_aminoacid.number
    end = x.last_aminoacid.number

file = file.replace(".pdb",f":{begin}:{end}_TF.pdb")
protein.write(file)

