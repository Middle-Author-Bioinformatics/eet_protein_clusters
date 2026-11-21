"""Deterministic whole-proteome extraction with explicit annotation provenance."""
from pathlib import Path
import re
from Bio import SeqIO
from Bio.SeqFeature import ExactPosition
from common import open_fasta, writer

MAP_FIELDS = ['protein_id', 'genome', 'original_id', 'original_header', 'annotation_source',
              'contig', 'start', 'end', 'strand', 'location', 'partial', 'translation_table', 'length_aa', 'gram']
NUCLEOTIDES = re.compile(r'[ACGTRYSWKMBDHVNU]+')


def discover(directory, input_type='auto'):
    suffixes = {'.fna': 'assembly', '.fa': 'assembly', '.fasta': 'assembly',
                '.gb': 'genbank', '.gbk': 'genbank', '.gbff': 'genbank'}
    if input_type == 'proteins': suffixes = {'.faa': 'proteins', '.fa': 'proteins', '.fasta': 'proteins'}
    files = []
    for path in sorted(directory.iterdir()):
        if not path.is_file(): continue
        bare = Path(path.name.removesuffix('.gz'))
        kind = suffixes.get(bare.suffix.lower())
        if not kind: continue
        if input_type not in ['auto', 'proteins'] and kind != input_type: continue
        name = bare.stem
        if any(name == existing[0] for existing in files):
            raise ValueError(f'Duplicate genome basename {name!r}; keep one annotation/input per genome.')
        files.append((name, path.resolve(), kind))
    if not files: raise ValueError(f'No {input_type} genome inputs found in {directory}; see GENOME_WORKFLOW.md')
    if input_type == 'auto' and any(p.name.endswith(('.faa', '.faa.gz')) for p in directory.iterdir()):
        raise ValueError('Protein FASTAs alongside genomes are ambiguous. Select --input-type assembly/genbank or --input-type proteins explicitly.')
    return files


def gram_assignments(files, default=None, metadata=None):
    import csv
    values = {}
    if metadata:
        with open(metadata) as f:
            reader = csv.DictReader(f, delimiter='\t')
            if not {'genome', 'gram'} <= set(reader.fieldnames or []):
                raise ValueError('Genome metadata requires tab-separated genome and gram columns')
            for row in reader:
                if row['genome'] in values: raise ValueError(f'Duplicate Gram metadata: {row["genome"]}')
                values[row['genome']] = row['gram']
        unknown = set(values)-{g for g, _, _ in files}
        if unknown: raise ValueError(f'Metadata names do not match input basenames: {sorted(unknown)}')
    result = {g: values.get(g, default) for g, _, _ in files}
    for genome, gram in result.items():
        if gram not in ['positive', 'negative']:
            raise ValueError(f'{genome}: supply --gram positive/negative or --genome-metadata; Gram class is never guessed.')
    return result


def records(path, kind, translation_table=11):
    """Yield (metadata, protein); assembly predictions never replace GenBank CDSs."""
    if kind == 'proteins':
        with open_fasta(path) as f:
            for rec in SeqIO.parse(f, 'fasta'):
                yield dict(original_id=rec.id, original_header=rec.description,
                           annotation_source='supplied_protein', partial='unknown'), str(rec.seq)
        return
    finder = None
    if kind == 'assembly':
        try: import pyrodigal
        except ImportError as exc: raise ValueError('Assembly gene prediction needs Pyrodigal: python -m pip install -r requirements-genomes.txt') from exc
        bins = pyrodigal.MetagenomicBins([b for b in pyrodigal.METAGENOMIC_BINS
                                        if b.training_info.translation_table == translation_table])
        if not bins: raise ValueError(f'No Pyrodigal metagenomic models for table {translation_table}; supported: 4,11')
        finder = pyrodigal.GeneFinder(meta=True, metagenomic_bins=bins, closed=False, mask=True)
    seen = set()
    with open_fasta(path) as f:
        for rec in SeqIO.parse(f, 'genbank' if kind=='genbank' else 'fasta'):
            if rec.id in seen: raise ValueError(f'{path}: duplicate contig ID {rec.id}')
            seen.add(rec.id)
            if not NUCLEOTIDES.fullmatch(str(rec.seq).upper()):
                raise ValueError(f'{path}/{rec.id}: invalid DNA; protein FASTAs require --input-type proteins')
            if finder:
                for i, gene in enumerate(finder.find_genes(str(rec.seq).upper()), 1):
                    yield dict(original_id=f'{rec.id}_{i}', original_header=f'{rec.id}_{i}',
                               annotation_source='pyrodigal_meta', contig=rec.id, start=gene.begin, end=gene.end,
                               strand='+' if gene.strand==1 else '-', location=f'{gene.begin}..{gene.end}',
                               partial=str(bool(gene.partial_begin or gene.partial_end)).lower(),
                               translation_table=translation_table), gene.translate(include_stop=False)
            else:
                for i, feature in enumerate(rec.features, 1):
                    if feature.type != 'CDS' or 'pseudo' in feature.qualifiers or 'pseudogene' in feature.qualifiers: continue
                    if feature.location is None: raise ValueError(f'{path}: CDS {i} lacks a location')
                    q = feature.qualifiers
                    original = q.get('protein_id', q.get('locus_tag', [f'{rec.id}_CDS{i}']))[0]
                    table = int(q.get('transl_table', [translation_table])[0])
                    source = 'genbank_translation'
                    if 'translation' in q: protein = q['translation'][0]
                    else:
                        if 'transl_except' in q: raise ValueError(f'{original}: translation exception without supplied translation; curate GenBank translation first')
                        dna = feature.extract(rec.seq)[int(q.get('codon_start', ['1'])[0])-1:]
                        dna = dna[:len(dna)//3*3]
                        protein = str(dna.translate(table=table, to_stop=False))
                        # Bacterial alternative initiators translate as M at a complete start.
                        from Bio.Data import CodonTable
                        start_pos = feature.location.parts[0].start if feature.location.strand != -1 else feature.location.parts[0].end
                        if isinstance(start_pos, ExactPosition) and str(dna[:3]) in CodonTable.unambiguous_dna_by_id[table].start_codons:
                            protein = 'M'+protein[1:]
                        source = 'genbank_CDS_translated'
                    partial = any(not isinstance(p, ExactPosition) for part in feature.location.parts for p in [part.start, part.end])
                    yield dict(original_id=original, original_header=original+' '+q.get('product', [''])[0],
                               annotation_source=source, contig=rec.id, start=int(feature.location.start)+1,
                               end=int(feature.location.end), strand={1:'+', -1:'-'}.get(feature.location.strand, '.'),
                               location=str(feature.location), partial=str(partial).lower(), translation_table=table), protein


def prepare(files, grams, output, translation_table=11):
    """Safe predictor IDs map to original genome/CDS IDs, including repeated IDs in different genomes."""
    output.mkdir(parents=True, exist_ok=True)
    proteomes = output/'proteomes'; proteomes.mkdir(exist_ok=True)
    for old in proteomes.glob('g*.faa'): old.unlink()
    sf, sw = writer(output/'protein_mapping.tsv', MAP_FIELDS)
    genomes = []
    try:
        with open(output/'proteins.faa', 'w') as combined:
            for gi, (genome, path, kind) in enumerate(files, 1):
                fasta = proteomes/f'g{gi:05d}.faa'; count = 0; seen = set()
                with open(fasta, 'w') as dest:
                    for count, (meta, seq) in enumerate(records(path, kind, translation_table), 1):
                        seq = str(seq).upper().rstrip('*')
                        if not seq or not re.fullmatch(r'[ACDEFGHIKLMNPQRSTVWYBXZJUO]+', seq):
                            raise ValueError(f'{genome}/{meta["original_id"]}: invalid protein or internal stop; curate annotation')
                        # Coordinate identity permits repeated GenBank protein IDs at different loci.
                        identity = (meta['original_id'], meta.get('contig', ''), meta.get('location', ''))
                        if identity in seen: raise ValueError(f'{genome}: duplicate protein/CDS {identity}')
                        seen.add(identity)
                        pid = f'g{gi:05d}p{count:08d}'
                        text = f'>{pid}\n{seq}\n'; dest.write(text); combined.write(text)
                        sw.writerow(dict(meta, protein_id=pid, genome=genome, gram=grams[genome], length_aa=len(seq)))
                if not count: raise ValueError(f'{genome}: no protein-coding sequences found; no other annotation mode was substituted')
                genomes.append(dict(genome=genome, input=str(path), input_type=kind, gram=grams[genome],
                                    proteins=count, fasta=str(fasta.resolve())))
    finally: sf.close()
    return genomes
