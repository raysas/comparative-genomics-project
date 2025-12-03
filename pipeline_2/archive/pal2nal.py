#!/usr/bin/env python3
import sys
import re

# -------------------------
# Help
# -------------------------
def show_help():
    print("""
Usage: pal2nal.py pep.aln nuc.fasta [nuc.fasta ...] [options] > output

Options:
    -output clustal|paml|fasta|codon
    -nogap
    -html
    -nostderr
    -blockonly
    -nomismatch
    -codontable universal|vmitochondria
""")
    sys.exit()

# -------------------------
# Parse arguments
# -------------------------
if len(sys.argv) < 3:
    show_help()

alnfile = None
nucfiles = []
outform = "clustal"
nogap = False
html = False
nostderr = False
blockonly = False
nomismatch = False
codontable = "universal"

i = 1
while i < len(sys.argv):
    arg = sys.argv[i]
    if arg == "-h":
        show_help()
    elif arg == "-output":
        i += 1
        outform = sys.argv[i]
        if outform not in ["clustal", "paml", "fasta", "codon"]:
            print("ERROR: valid output format: clustal, paml, fasta, or codon", file=sys.stderr)
            sys.exit(1)
    elif arg == "-nogap":
        nogap = True
    elif arg == "-html":
        html = True
    elif arg == "-nostderr":
        nostderr = True
    elif arg == "-blockonly":
        blockonly = True
    elif arg == "-nomismatch":
        nomismatch = True
    elif arg == "-codontable":
        i += 1
        codontable = sys.argv[i]
        if codontable not in ["universal", "vmitochondria"]:
            print("ERROR: valid codon table: universal or vmitochondria", file=sys.stderr)
            sys.exit(1)
    elif not alnfile:
        alnfile = arg
    else:
        nucfiles.append(arg)
    i += 1

# -------------------------
# Read nucleotide sequences
# -------------------------
id2nucseq = {}
for file in nucfiles:
    with open(file, "r") as f:
        data = f.read()
    data = re.sub(r'\r\n|\r|\n', '\n', data)
    tmpid = None
    for line in data.split('\n'):
        if not line.strip() or line.startswith('#'):
            continue
        if line.startswith('>'):
            tmpid = line[1:].split()[0]
            id2nucseq[tmpid] = ""
        elif tmpid:
            seq = re.sub(r'[^a-zA-Z]', '', line)
            id2nucseq[tmpid] += seq

# -------------------------
# Read protein alignment
# -------------------------
with open(alnfile, "r") as f:
    aln_data = f.read()
aln_data = re.sub(r'\r\n|\r|\n', '\n', aln_data)
lines = aln_data.split('\n')

# Detect format
inalntype = "clustal"
for line in lines:
    if line.startswith('CLUSTAL'):
        inalntype = "clustal"
        break
    elif line.startswith('>'):
        inalntype = "fasta"
        break
    elif line.startswith('Gblocks'):
        inalntype = "gblocks"
        break

# Parse sequences
id2aaaln = {}
aaid = []

if inalntype == "fasta":
    tmpid = None
    for line in lines:
        if line.startswith('>'):
            tmpid = line[1:].split()[0]
            aaid.append(tmpid)
            id2aaaln[tmpid] = ""
        else:
            seq = re.sub(r'\s+', '', line).upper()
            if tmpid:
                id2aaaln[tmpid] += seq

elif inalntype == "clustal":
    for line in lines:
        if not line.strip() or line.startswith('CLUSTAL'):
            continue
        if re.match(r'^\S+', line):
            parts = line.split()
            seqid = parts[0]
            seqpart = parts[1].upper()
            if seqid not in id2aaaln:
                id2aaaln[seqid] = ""
                aaid.append(seqid)
            id2aaaln[seqid] += seqpart

elif inalntype == "gblocks":
    # Simple Gblocks parsing (like fasta)
    tmpid = None
    for line in lines:
        if line.startswith('>'):
            tmpid = line[1:].split()[0]
            aaid.append(tmpid)
            id2aaaln[tmpid] = ""
        else:
            seq = re.sub(r'\s+', '', line).upper()
            if tmpid:
                id2aaaln[tmpid] += seq

# -------------------------
# Build codon alignment
# -------------------------
id2codonseq = {}
for aid in aaid:
    aa_seq = id2aaaln[aid]
    nuc_seq = id2nucseq.get(aid)
    if not nuc_seq:
        if not nostderr:
            print(f"Warning: no nucleotide sequence for {aid}", file=sys.stderr)
        id2codonseq[aid] = ""
        continue

    codons = [nuc_seq[i*3:(i+1)*3] for i in range(len(aa_seq))]
    codon_aln = ""
    for aa, codon in zip(aa_seq, codons):
        if aa == "-":
            codon_aln += "---"
        else:
            if len(codon) != 3:
                codon = "NNN"
            codon_aln += codon
    id2codonseq[aid] = codon_aln

# -------------------------
# Remove gaps if requested
# -------------------------
if nogap:
    for aid in id2codonseq:
        codseq = id2codonseq[aid]
        codseq = "".join([codseq[i:i+3] for i in range(0, len(codseq), 3) if id2aaaln[aid][i//3] != "-"])
        id2codonseq[aid] = codseq

# -------------------------
# Output
# -------------------------
def write_fasta(seqdict):
    for seqid, seq in seqdict.items():
        print(f">{seqid}")
        print(seq)

def write_clustal(seqdict):
    print("CLUSTAL W alignment\n")
    keys = list(seqdict.keys())
    length = len(next(iter(seqdict.values())))
    block_size = 60
    for start in range(0, length, block_size):
        for k in keys:
            print(f"{k:<15} {seqdict[k][start:start+block_size]}")
        print("")

if outform == "fasta":
    write_fasta(id2codonseq)
elif outform == "clustal":
    write_clustal(id2codonseq)
elif outform == "codon":
    write_fasta(id2codonseq)
elif outform == "paml":
    # simple PAML-like format
    nseqs = len(id2codonseq)
    seqlen = len(next(iter(id2codonseq.values())))
    print(f" {nseqs} {seqlen}")
    for k, v in id2codonseq.items():
        print(f"{k} {v}")