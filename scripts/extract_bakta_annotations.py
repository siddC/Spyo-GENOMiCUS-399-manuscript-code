#!/usr/bin/env python3
"""
Extract AMR, virulence (VFDB), IS-element and FCT/pilin-locus gene content from
Bakta JSON outputs for Spyo-GENOMiCUS-399.

Usage:
    python3 extract_bakta_annotations.py /path/to/bakta_dir  out_prefix

Expects one subdirectory per genome, each containing <genome_id>.json
(the Bakta layout).  Writes two small TSVs:

    <out_prefix>_expert_long.tsv   one row per expert-system hit
    <out_prefix>_fct_long.tsv      one row per FCT/pilin-locus gene
"""
import sys, os, json, glob, re, csv

# Acquired-AMR nomenclature is parenthesised -- tet(O), erm(B), mef(A), msr(D).
# This distinguishes them from tetR regulators and msrA/B/C methionine
# sulfoxide reductases, which are NOT resistance determinants.
AMR_UNPARENTHESISED = {
    'catq','cat','dfrg','dfra','sat4','aph(3\')-iii','ant(6)-ia','aad6',
    'vana','vanb','blaz','fosa','fosb','optra','poxta','cfr',
}

FCT = {
    # FCT-region pilin/sortase genes plus the mga-regulon surface-protein genes.
    # Reported as gene CONTENT only -- this is not formal FCT typing.
    # 'lepA' is deliberately excluded: it collides with translation elongation factor 4.
    'cpa','fcta','fctb','srtc1','srtc2','srtc','sipa','tee',
    'prtf1','sfbi','prtf2','sfbii','nra','rofa','msmr','mga','emm','enn','mrp',
}

def norm(s):
    return re.sub(r'[^a-z0-9()\'-]', '', str(s or '').lower())

def is_amr_symbol(g):
    g = norm(g)
    if re.match(r"^(tet|erm|mef|msr|lnu|lsa|vga|vat|aac|aph|ant|aad|cfr)\(", g):
        return True
    return g in AMR_UNPARENTHESISED

def main(bakta_dir, prefix):
    paths = sorted(glob.glob(os.path.join(bakta_dir, '*', '*.json')))
    if not paths:
        sys.exit(f'No <genome>/<genome>.json found under {bakta_dir}')
    ew = csv.writer(open(f'{prefix}_expert_long.tsv', 'w', newline=''), delimiter='\t')
    fw = csv.writer(open(f'{prefix}_fct_long.tsv',    'w', newline=''), delimiter='\t')
    ew.writerow(['genome_id','category','locus','gene','product','identity',
                 'query_cov','subject_cov','method','ref_id','contig','start','stop','strand'])
    fw.writerow(['genome_id','locus','gene','product','contig','start','stop','strand'])

    n_amr = n_gen = 0
    for p in paths:
        gid = os.path.splitext(os.path.basename(p))[0]
        try:
            with open(p) as fh:
                d = json.load(fh)
        except Exception as e:
            print(f'  ! skipped {gid}: {e}', file=sys.stderr); continue
        n_gen += 1
        for ft in d.get('features', []):
            locus  = ft.get('locus') or ''
            contig = ft.get('sequence') or ''
            st, sp, sd = ft.get('start'), ft.get('stop'), ft.get('strand')

            for e in (ft.get('expert') or []):
                # AMRFinderPlus entries carry type=='amrfinder' and NO 'source' key.
                src = e.get('source')
                cat = 'amrfinder' if e.get('type') == 'amrfinder' else (src or 'unknown')
                if cat == 'amrfinder':
                    n_amr += 1
                elif cat not in ('VFDB', 'IS', 'BlastRules'):
                    continue
                ref = next((x.split(':',1)[1] for x in (e.get('db_xrefs') or [])), e.get('id',''))
                ew.writerow([gid, cat, locus, e.get('gene') or ft.get('gene') or '',
                             e.get('product') or '', e.get('identity',''),
                             e.get('query_cov',''), e.get('subject_cov',''),
                             e.get('method',''), ref, contig, st, sp, sd])

            names = {norm(ft.get('gene'))} | {norm(g) for g in (ft.get('genes') or [])}
            if names & FCT:
                fw.writerow([gid, locus, ft.get('gene') or '', ft.get('product') or '',
                             contig, st, sp, sd])

    print(f'Parsed {n_gen} genomes; {n_amr} AMRFinderPlus hits total.')
    print(f'Wrote {prefix}_expert_long.tsv and {prefix}_fct_long.tsv')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
