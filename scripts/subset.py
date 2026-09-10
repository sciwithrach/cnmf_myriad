# import packages
import scanpy as sc
import anndata as ad
import argparse
from pathlib import Path

def build_parser():
	
	p = argparse.ArgumentParser(
        	description=__doc__,
        	formatter_class=argparse.RawDescriptionHelpFormatter,
	)

	p.add_argument(
		"--counts", 
		metavar="FILE",
                help="AnnData file (.h5ad)"
	)
	p.add_argument(
		"--metadata", 
                help="Metadata column in obs e.g. celltype"
	)
	p.add_argument(
		"--subset", 
        	help="Subset of metadata to keep e.g. HC"
		)
	p.add_argument(
		'--output-dir',
                help="top-level directory for saving outputs"
	)
    
	return p

def main():
	
	# setup arguments
	args = build_parser().parse_args()
	
	# setup outpath
	outpath = Path(args.output_dir)
	outpath.mkdir(parents=True, exist_ok=True)
	# testing
	print(outpath)
	
	# read adata
	adata = ad.read_h5ad(args.counts)

	# subset by obs
	adata = adata[adata.obs[args.metadata] == args.subset]

	# write adata
	adata.write(f'{outpath}/adata_{args.metadata}_{args.subset}.h5ad')

if __name__ == "__main__":
    main()

