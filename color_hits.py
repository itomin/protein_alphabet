# from pymol import cmd
# from pymol import stored
#
#
# def color_hit(pdb_id, start, end):
#     print("fetch")
#     cmd.fetch(pdb_id, async=0)
#     print("(resi %s-%s)" % (start, end))
#     cmd.select("m1", "(resi %s-%s)" % (start, end))
#     # cmd.select("m1", "(resi 1-10)")
#     cmd.color("purple", "m1")
#
#
#
# cmd.extend("color_hit",color_hit)

# run "C:\Users\tomin\dev\suffix_alphabet\color_hits.py"
# load C:\Users\tomin\dev\data\pdb_str\pdb\12as.ent


import pymol
pymol.finish_launching()

def gen_pymol_session(pdb1, pdb2, pos1, pos2):
    pymol.cmd.reinitialize()
    pymol.cmd.fetch(pdb1, async_=0)
    pymol.cmd.color("palegreen", pdb1)
    pymol.cmd.fetch(pdb2, async_=0)
    pymol.cmd.color("lightblue", pdb2)
    pymol.cmd.select("frag1_sele", pdb1 + " and resi %s-%s" % pos1)
    pymol.cmd.select("frag2_sele", pdb2 + " and resi %s-%s" % pos2)
    pymol.cmd.color("green", "frag1_sele")
    pymol.cmd.color("blue", "frag2_sele")
    pymol.cmd.super("frag1_sele", "frag2_sele")
    pymol.cmd.copy_to("frag1_sele", "frag1")
    pymol.cmd.copy_to("frag2_sele", "frag2")
    pymol.cmd.hide("everything")
    pymol.cmd.show("cartoon")


def check_selection_resolved(pdb, sele):
    pymol.cmd.reinitialize()
    pymol.cmd.fetch(pdb, async_=0)
    resis_resolved = pymol.cmd.iterate_state(-1, sele + " and name ca", 'pass')
    resis = pymol.cmd.iterate(sele + " and name ca", 'pass')
    print(resis_resolved)
    print(resis)
    if resis_resolved == resis:
        return True
    else:
        return False



print(check_selection_resolved("13gs", "13gs"))
# print(check_selection_resolved("12as", "12as and resi 6-20"))

# gen_pymol_session("1aro", "2zsg", (676,690), (287,301) )
# pymol.cmd.save("sessions/chem/1aro_2zsg.pse")
#
# gen_pymol_session("1at9", "5g0x", (92,106), (120,134))
# pymol.cmd.save("sessions/chem/1at9_5g0x.pse")
#
# gen_pymol_session("1e8c", "5t5q", (137,151), (9,23))
# pymol.cmd.save("sessions/chem/1e8c_5t5q.pse")
#
# gen_pymol_session("1aor", "3opy", (102,116), (138,152))
# pymol.cmd.save("sessions/chem/1aor_3opy.pse")
#
# gen_pymol_session("1ddq", "3goh", (848,858), (64,74))
# pymol.cmd.save("sessions/chem/1ddq_3goh.pse")
#
# gen_pymol_session("1ci9", "5kjp", (184,198), (187,201))
# pymol.cmd.save("sessions/chem/1ci9_5kjp.pse")
#
# gen_pymol_session("1bg5", "2zxe", (237,251), (157,171))
# pymol.cmd.save("sessions/chem/1bg5_2zxe.pse")
#
# gen_pymol_session("1chu", "2l35", (10,24), (37,51))
# pymol.cmd.save("sessions/chem/1chu_2l35.pse")
#
# gen_pymol_session("1ddq", "1hi9", (576,590), (146,160))
# pymol.cmd.save("sessions/chem/1ddq_1hi9.pse")
#
# gen_pymol_session("1bkb", "2q2b", (82,92), (221,231))
# pymol.cmd.save("sessions/chem/1bkb_2q2b.pse")
#
# gen_pymol_session("1ci4", "1ehk", (19,30), (195,206))
# pymol.cmd.save("sessions/chem/1ci4_1ehk.pse")
#
# gen_pymol_session("1df4", "2cmr", (3,33), (7,37))
# pymol.cmd.save("sessions/chem/1df4_2cmr.pse")
#
# gen_pymol_session("1at9", "1bha", (37,70), (35,68))
# pymol.cmd.save("sessions/chem/1at9_1bha.pse")

# gen_pymol_session("1a2p", "3da7", (7,66), (46,105))
# pymol.cmd.save("sessions/chem/1a2p_3da7.pse")
