# import boto
# from boto.s3.key import Key
# from boto.s3.connection import S3Connection
import pandas as pd
from Bio import SeqIO
# from Bio.PDB import PDBList
# from functools import reduce
import seaborn as sns
import matplotlib.pyplot as plt

%matplotlib inline

keyId ="AKIAJ4XU3KONWZ7BCIUA"
sKeyId="g7rUiUhA4Weu5V5A9hiok7Ftj+FQuTemwEXBZBWJ"
REGION_HOST =  's3.eu-central-1.amazonaws.com'
bucketName="cathfasta"

# dbl = Bio.PDB.PDBList()
# ff = dbl.retrieve_pdb_file("2a50


def extract_interval(x):
    _from, to_chain = x.split("-")
    _to, chain = to_chain.split(":")
    print([_from, _to])
    return [_from, _to]


# read cs219
def readCS219(file):
    with open("../data/pdb_str/cs219/%s.as" % file, "r",  encoding="cp437") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            cs219 = record.seq
    return cs219



# read fasta
def readPDB(file):
    for record in SeqIO.parse("../data/pdb_str/fasta/%s.fasta" % file, "fasta"):
        pdb = record.seq
    return pdb


# conn = S3Connection(keyId,sKeyId, REGION_HOST)

#---------------------------------------------------------------
#
# retrieve all fasta files
#
#---------------------------------------------------------------
# bucket = conn.get_bucket("cathfasta")
# for key in bucket:
#     if key.name.endswith(".fasta"):
#         key.get_contents_to_filename( '../data/pdb_str/fasta/' + key.name)


#---------------------------------------------------------------
#
# retrieve all cs219 files
#
#---------------------------------------------------------------
# bucket = conn.get_bucket("cathcs219")
# for key in bucket:
#     key.get_contents_to_filename( '../data/pdb_str/cs219/' + key.name)


pdb_intervals = pd.read_csv("data/cath-b-newest-all", sep=" ", names =["pdb", "version", "cath", "range"])
pdb_intervals = pdb_intervals[["pdb", "range"]]

df = pd.read_csv("../data/pdb_str/part-r-00000-be9d15c5-c7dd-48bd-bbc1-27879ce81b6a.csv",  \
                sep=",", escapechar="\\", header=None,  \
                names=["hit", "pdb", "c", "a", "t", "h", "N"])


# 1gpw
# 1vh7

# df = df[df.N > 10]
# df = df[df.pdb == "2a50:C00"]




df.loc[:,"pdb"] = df["pdb"].apply(lambda x: "".join(x.split(":")))
df.loc[:,"pdb_seq"]  = df["pdb"].apply(lambda x: readPDB(x))
df.loc[:,"cs219"]  = df["pdb"].apply(lambda x: readCS219(x))

df.loc[:,"pdb_seq_len"] = df["pdb_seq"].apply(lambda x: len(x))
df.loc[:,"cs219_len"] = df["cs219"].apply(lambda x: len(x))

df.loc[:, "start"] = df[["cs219", "hit"]].apply(lambda x: x.cs219.find(x.hit) + 1, axis=1)
df.loc[:, "end"] = df[["start", "hit"]].apply(lambda x: x.start + len(x.hit) - 1, axis=1)


df = pd.merge(df, pdb_intervals, how="inner", on="pdb")


#TODO offset bestimmen, pdb chain betrachten
#TODO a3m in msa konvertieren
# index offset aus pdb laden, 1krl_A00 beginnt nicht bei 1 sondern 101
# pymol script schreiben


clearEmpty = lambda x: [ i for i in x if i != '']
def parse_interval(range):
    segs = [ seg.split("-")  for seg in range.split(",")]
    segs = [ clearEmpty(s) for s in segs ]
    return [ [ int(s[0]), int(s[1].split(":")[0]) ] for s in segs]



def merge_range(elem):
    return ",".join([ ":".join( [str(o) for o in pair] ) for pair in elem])


def scale_hit(pdb_seq, start, end):
    cs219_interval = []
    # diffs_start = 0
    num_char_hit = (end - start + 1)
    num_char_viewed_hit = 0
    num_char_viewed_seq = 0

    for I in pdb_seq:
        start_scaled = max((I[0] + start - 1) - num_char_viewed_seq, I[0])
        # end_scaled = min(I[0] + start + num_char_hit - num_char_viewed_seq, I[1])
        end_scaled = min(start_scaled+num_char_hit-1, I[1])
        # print("num_char_viewed_seq %s" % num_char_viewed_seq)
        # print("num_char_viewed_hit %s" % num_char_viewed_hit)
        # print("num_char_hit %s" % num_char_hit)
        # print("start_scaled %s " % start_scaled)
        # print("end_scaled %s " % end_scaled)
        if start_scaled <= I[1] and num_char_hit > 0:
            cs219_interval = cs219_interval + [[start_scaled, end_scaled]]
            # diffs_end = diffs_end + end_scaled - start_scaled + 1
            num_char_hit = num_char_hit - (end_scaled-start_scaled+1)
            num_char_viewed_hit = end_scaled - start_scaled + 1
            num_char_viewed_seq = num_char_viewed_seq + (start_scaled-start)+num_char_viewed_hit # + end_scaled - start
            # start = 1
        else:
            num_char_viewed_seq = num_char_viewed_seq +(I[1]-I[0])+1

    return cs219_interval

# pdb_seq = [[3,6],[10,20], [25, 30]]
# scale_hit(pdb_seq, 4, 7) # 6,6 - 10,12
# scale_hit(pdb_seq, 1, 7) # 3,6 - 10,12
# scale_hit(pdb_seq, 4, 4) # 6,6
# scale_hit(pdb_seq, 3, 4) # 5,6
# scale_hit(pdb_seq, 5, 7) # 10,12
#
# scale_hit(pdb_seq, 5, 16) # 10,20 - 25
# scale_hit(pdb_seq, 7, 16)



df.loc[:, "range_scaled"] = df[["range", "start", "end"]].apply(lambda x: merge_range(scale_hit(parse_interval(x.range), x.start, x.end)), axis=1)
df.head(2)

hits = df.groupby(["hit"]).size().reset_index(name="Number_of_Seq")
hits["index"] = hits.index


df = pd.merge(df, hits, how="inner", on="hit")

import csv
df.to_csv("cs219_found_hits.csv", index=False, quoting=csv.QUOTE_NONNUMERIC, escapechar="\\")


# which sequences are matched that are not in the same class, in the same architecture?

hits = df[["hit", "c", "N"]].drop_duplicates().groupby(["hit", "N"]).size().reset_index(name="freq")
hits = hits[hits.freq > 1][["hit"]]



pd.merge(hits, df, how="inner", on="hit")[["hit", "N", "pdb", "start", "end", "c", "a", "t", "pdb_seq_len", "cs219_len", "range", "range_scaled" ]]#.to_csv('../data/pdb_str/pymol/sample.csv')


readCS219("2c9jA00")
readPDB("2fl1A00")[16:34]
readPDB("2c9jA00")[9:27]
