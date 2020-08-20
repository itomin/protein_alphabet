import pymol
import pandas as pd
import logging
import math
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.INFO)
logger = logging.getLogger('logger')


def super_frag(file, x,y, x_frag, y_frag):
        x_orig = x
        y_orig = y
        x = x[:4] + "_" + x[4]
        y = y[:4] + "_" + y[4]

        def get_first_resi(sele):
            myspace = {'first_resi': [] }
            pymol.cmd.iterate('first ' + sele, 'first_resi.append(resi)', space=myspace)
            return(int(myspace['first_resi'][0]))

        logger.debug("PyMol: Reinitialize")
        pymol.cmd.reinitialize()
        logger.debug("PyMol: Fetch first pdb")
        pymol.cmd.fetch(x)
        logger.debug("PyMol: Fetch second pdb")
        pymol.cmd.fetch(y)
        pymol.cmd.remove("solvent")
        pymol.cmd.color("lightorange", x)
        pymol.cmd.color("palegreen", y)

        seq_fetched_1 = "".join(pymol.cmd.get_fastastr(x).splitlines()[1:])
        seq_fetched_2 = "".join(pymol.cmd.get_fastastr(y).splitlines()[1:])

        hit_start_1 = seq_fetched_1.find(x_frag) + get_first_resi(x)
        hit_end_1 =  hit_start_1 + len(x_frag) - 1
        hit_start_2 = seq_fetched_2.find(y_frag) + get_first_resi(y)
        hit_end_2 =  hit_start_2 + len(y_frag) - 1

        logger.debug("Hit in %s: %d-%d, Hit in %s: %d-%d" % (x, hit_start_1, hit_end_1,
                                                             y, hit_start_2, hit_end_2))
        def check_selection_resolved(sele):
            resis_resolved= pymol.cmd.iterate_state(-1, sele + " and name ca", 'pass')
            resis = pymol.cmd.iterate(sele + " and name ca", 'pass')

            if resis_resolved == resis:
                return True
            else:
                return False

        frag_start_sele = "%s_fragment_sele" % (x)
        frag_start = "%s_fragment" % (x)
        pymol.cmd.select(frag_start_sele, "%s and resi %d-%d" % (x, hit_start_1, hit_end_1))
        frag_end_sele = "%s_fragment_sele" % (y)
        frag_end = "%s_fragment" % (y)
        pymol.cmd.select(frag_end_sele, "%s and resi %d-%d" % (y, hit_start_2, hit_end_2))

        logger.debug(frag_end_sele + ": %d-%d" % (hit_start_2, hit_end_2))
        logger.debug(frag_start_sele + ": %d-%d" % (hit_start_1, hit_end_1))

        rmsd = pymol.cmd.super(frag_start_sele, frag_end_sele)
        pymol.cmd.color("orange", frag_start_sele)
        pymol.cmd.color("green", frag_end_sele)
        pymol.cmd.copy_to(frag_start, frag_start_sele)
        pymol.cmd.copy_to(frag_end, frag_end_sele)

        pymol.cmd.group("fragment", " or ".join([x,
                                                 y,
                                                frag_start,
                                                frag_end,
                                                frag_start_sele,
                                                frag_end_sele]))

        pymol.cmd.save("sessions/%s/%s_%s.pse" % (file,x_orig,y_orig))

#----------------------------------------------------------------------------
# load all aa alphabets
#----------------------------------------------------------------------------
aa = pd.read_csv("hit_aa_sample.csv")
aa = aa[aa.hit_length >= 10]
aa.head()


for idx, row in aa.iterrows():
    try:
        super_frag("aa", row.x, row.y, row.hit, row.hit)
    except:
        print("failed %s_%s" % (row.x, row.y))

#----------------------------------------------------------------------------
# load all chem alphabets
#----------------------------------------------------------------------------
chem = pd.read_csv("chem_hits_sample.csv")
chem.head()

super_frag("chem", "1e5mA01", "1e8cA02", "IGVLIGTGIGGL", "LLLALATLLALG")

for idx, row in chem.iterrows():
    try:
        super_frag("chem", row.x, row.y, row.x_pdb_fragment, row.y_pdb_fragment)
    except:
        print("failed %s_%s" % (row.x, row.y))


#----------------------------------------------------------------------------
# load all chem alphabets
#----------------------------------------------------------------------------
chem = pd.read_csv("chem_hits_sample.csv")
chem.head()

for idx, row in chem.iterrows():
    try:
        super_frag("chem", row.x, row.y, row.x_pdb_fragment, row.y_pdb_fragment)
    except:
        print("failed %s_%s" % (row.x, row.y))


#----------------------------------------------------------------------------
# Load all pb alphabets
#-----------------------------------------------------------------------------
pb = pd.read_csv("pb_hits_sample.csv")
pb.head(2)


for idx, row in pb.iterrows():
    try:
        super_frag("pb", row.x, row.y, row.x_pdb_fragment, row.y_pdb_fragment)
    except:
        print("failed %s_%s" % (row.x, row.y))



#----------------------------------------------------------------------------
# Load all m32k2 alphabets
#-----------------------------------------------------------------------------
m32k25 = pd.read_csv("m32k2__v2_hits_sample.csv")
m32k25.head(2)


for idx, row in m32k25.iterrows():
    try:
        super_frag("m32k25_v2", row.x, row.y, row.x_pdb_fragment, row.y_pdb_fragment)
    except:
        print("failed %s_%s" % (row.x, row.y))
