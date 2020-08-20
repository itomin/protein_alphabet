import pymol
import pandas as pd
import logging
import math
logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.ERROR)
logger = logging.getLogger('logger')

filter_frag = pd.read_csv("hit_aa_sample.csv")
filter_frag = filter_frag[filter_frag.hit_length >= 10]

# fragments = pd.read_csv("chem_hits_sample.csv")
fragments = pd.read_csv("m32k2_hits_sample.csv")
# fragments = fragments[fragments.x_pdb_fragment != fragments.y_pdb_fragment]
# fragments = fragments[fragments.hit_length > 15]
fragments = fragments[fragments.m32k25_fragment.apply(lambda x: len(set(x)) > 5)]
fragments = fragments.sample(150,replace=False )
fragments = pd.merge(fragments, filter_frag, on=["x","y"], how="left")

fragments = fragments[fragments.hit.apply(lambda x: pd.isnull(x))]

fragments = fragments[~ fragments.x.str.startswith("153l")]
fragments.shape

def gen_pymol_session(frag, aa=False):
        start_node = frag[0][:4] + "_" + frag[0][4]
        end_node = frag[1][:4] + "_" + frag[1][4]
        start_hit_seq = frag[4]
        end_hit_seq = frag[5]
        def get_first_resi(sele):
            myspace = {'first_resi': [] }
            pymol.cmd.iterate('first ' + sele, 'first_resi.append(resi)', space=myspace)
            return(int(myspace['first_resi'][0]))

        logger.debug("PyMol: Reinitialize")
        pymol.cmd.reinitialize()
        logger.debug("PyMol: Fetch first pdb")
        pymol.cmd.fetch(start_node)
        logger.debug("PyMol: Fetch second pdb")
        pymol.cmd.fetch(end_node)
        pymol.cmd.remove("solvent")
        pymol.cmd.color("lightorange", start_node)
        pymol.cmd.color("palegreen", start_node)


        seq_fetched_1 = "".join(pymol.cmd.get_fastastr(start_node).splitlines()[1:])
        seq_fetched_2 = "".join(pymol.cmd.get_fastastr(end_node).splitlines()[1:])

        hit_start_1 = seq_fetched_1.find(start_hit_seq) +1 + get_first_resi(start_node)
        hit_end_1 =  hit_start_1 + len(start_hit_seq)
        hit_start_2 = seq_fetched_2.find(end_hit_seq) +1 + get_first_resi(end_node)
        hit_end_2 =  hit_start_2 + len(end_hit_seq)

        logger.debug("Hit in %s: %d-%d, Hit in %s: %d-%d" % (start_node, hit_start_1,
                                                            hit_end_1,end_node,
                                                            hit_start_2, hit_end_2))
        def check_selection_resolved(sele):
            resis_resolved= pymol.cmd.iterate_state(-1, sele + " and name ca", 'pass')
            resis = pymol.cmd.iterate(sele + " and name ca", 'pass')

            if resis_resolved == resis:
                return True
            else:
                return False

        frag_start_sele = "%s_fragment_sele" % (start_node,)
        frag_start = "%s_fragment" % (start_node,)
        pymol.cmd.select(frag_start_sele, "%s and resi %d-%d" % (start_node, hit_start_1, hit_end_1))
        frag_end_sele = "%s_fragment_sele" % (end_node,)
        frag_end = "%s_fragment" % (end_node,)

        pymol.cmd.select(frag_end_sele, "%s and resi %d-%d" % (end_node, hit_start_2, hit_end_2))

        logger.debug(frag_end_sele + ": %d-%d" % (hit_start_2, hit_end_2))
        logger.debug(frag_start_sele + ": %d-%d" % (hit_start_1, hit_end_1))

        rmsd = pymol.cmd.super(frag_start_sele, frag_end_sele)
        pymol.cmd.color("orange", frag_start_sele)
        pymol.cmd.color("green", frag_end_sele)
        pymol.cmd.copy_to(frag_start, frag_start_sele)
        pymol.cmd.copy_to(frag_end, frag_end_sele)

        pymol.cmd.group("fragment", " or ".join([start_node,
                                                        end_node,
                                                        frag_start,
                                                        frag_end,
                                                        frag_start_sele,
                                                        frag_end_sele]))

        pymol.cmd.save("sessions/m32k25/%s_%s.pse" % (start_node, end_node))

fragments.head()


for frag in fragments.iterrows():
    try:
        gen_pymol_session(frag[1])
    except:
        print("Failed to generate PyMol Session")
        continue
