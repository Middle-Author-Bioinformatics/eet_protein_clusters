#!/usr/bin/env python3
"""Generate artificial assemblies and GA-bearing models; no biological validation claims."""
import argparse
from pathlib import Path
import random
import json
from Bio.Data import CodonTable
import pyhmmer as ph


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', type=Path, required=True)
    args=p.parse_args()
    if args.output.exists() and any(args.output.iterdir()): p.error('Choose an empty output directory')
    (args.output/'genomes').mkdir(parents=True,exist_ok=True)
    rng = random.Random(42); alphabet=ph.easel.Alphabet.amino()
    # Biased but diverse coding sequences yield reproducible bacterial ORF calls.
    codons = ['GCT','GAA','AAA','GGT','GTT','GAC','TCT','GCC','CTG','ACG','TTC','CGT','CAT','TGT','ATC','AAC','CAG','TAT','TGG','CCT']
    table=CodonTable.unambiguous_dna_by_id[11]
    coding = {name: ''.join(rng.choices(codons,k=180)) for name in ['A','B']}
    aa = {name: ''.join(table.forward_table[s[i:i+3]] for i in range(0,len(s),3)) for name,s in coding.items()}
    with open(args.output/'models.hmm', 'wb') as f:
        for name, sequence in aa.items():
            msa=ph.easel.TextMSA(name='Demo_'+name, sequences=[ph.easel.TextSequence(name='template',sequence=sequence)])
            hmm, _, _ = ph.plan7.Builder(alphabet).build_msa(msa.digitize(alphabet), ph.plan7.Background(alphabet))
            hmm.accession='DEMO_'+name+'.1';hmm.cutoffs.gathering=(20.,20.);hmm.write(f)
    motif_dna='TGTGCTGCTTGTCAT'+'GGT'*3+'TGTGCTGCTTGTCAT'+'GGT'*3+'TCTGCTTCTCCTACT'
    for gi in [1,2]:
        blocks = []
        for i in range(12):
            orf='ATG'+coding['A' if i%3 else 'B']+(motif_dna if i%4==0 else '')+'TAA'
            blocks.append('TAA'+'A'*15+'AGGAGG'+'A'*7+orf+'C'*30)
        (args.output/'genomes'/f'demo_{gi}.fna').write_text(f'>synthetic_{gi}\n'+''.join(blocks)+'\n')
    (args.output/'MODEL_VERSION.txt').write_text('Artificial demo A/B models; seed 42; GA sequence/domain 20 bits. Not Pfam or validated EET models.\n')
    print(json.dumps(dict(genomes=str(args.output/'genomes'),hmm=str(args.output/'models.hmm'),localization='not supplied; explicitly skip or install licensed predictors'),indent=2))


if __name__=='__main__': main()
