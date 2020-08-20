import os

HOME = "/home/tomin/"
atoms = "atoms"

# for f in *; do mv "$f" "${f/%/.pdb}"; done

for F in os.listdir(HOME + "atoms"):
    os.system("g_sa_encode -s atoms/%s -f atoms/%s -strlf m32k25/%s  -rmsdlf tmp.rmsd -xpmlf tmp.lf -alphabet M32K25 -log log.log" % (F,F,F))
    os.system("rm *log*")
    os.system("rm *rmsd*")
