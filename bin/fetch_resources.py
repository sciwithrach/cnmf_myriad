#!/usr/bin/env python3
"""Download the reference files the GEP analysis needs, once, on a machine with internet (e.g. a login node).

    python bin/fetch_resources.py [--outdir assets/resources]

Writes MOUSE_GO_bp_no_GO_iea_symbol.gmt (Bader lab GO biological process gene sets) and collectri_mouse.parquet
(CollecTRI transcription factor network, translated to mouse by decoupler). Compute nodes may not have internet
access, so the analysis reads these files instead of downloading anything.
"""
import argparse
import urllib.request
from pathlib import Path

GMT_URL = 'https://download.baderlab.org/EM_Genesets/current_release/Mouse/symbol/GO/MOUSE_GO_bp_no_GO_iea_symbol.gmt'


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--outdir', default='assets/resources')
    args = p.parse_args()
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    gmt = out / GMT_URL.rsplit('/', 1)[1]
    if gmt.exists():
        print(f'{gmt} already exists, skipping')
    else:
        print(f'Downloading {GMT_URL}')
        urllib.request.urlretrieve(GMT_URL, gmt)

    tf = out / 'collectri_mouse.parquet'
    if tf.exists():
        print(f'{tf} already exists, skipping')
    else:
        import decoupler as dc
        dc.op.collectri(organism='mouse', verbose=True).to_parquet(tf)
        print(f'Wrote {tf}')


if __name__ == '__main__':
    main()
