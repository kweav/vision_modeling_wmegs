#!/usr/bin/env python3

import argparse
import numpy as np
import pandas as pd


def main():
	parser = generate_parser()
	args = parser.parse_args()
	outdf = concatenate_together(args)
	outdf.to_csv(args.outfile, sep ='\t', mode='a')


def generate_parser():
	parser = argparse.ArgumentParser(description = "Concatenate output files from Joint EP modeling to make a long dataframe with Condition, Genome, Features, CT, Rep, and R2")
	parser.add_argument('--standard-control-both', dest='scb', type=str, action='store', required=True, help='output file from standard cre refinement, no shuffling, and both features used in prediction')
	parser.add_argument('--standard-shufflecre-both', dest='sscb', type=str, action='store', required=True, help='output file from standard cre refinement, cre shuffling, and both features used in prediction')
	parser.add_argument('--standard-shuffletss-both', dest='sstb', type=str, action='store', required=True, help='output file from standard cre refinement, tss shuffling, and both features used in prediction')
	parser.add_argument('--standard-control-creonly', dest='scco', type=str, action='store', required=True, help='output file from standard cre refinemenet, no shuffling, and only cre features used in prediction')
	parser.add_argument('--standard-control-tssonly', dest='scto', type=str, action='store', required=True, help='output file from standard cre refinement, no shuffling, and only tss features used in prediction')
	parser.add_argument('--standard-shufflecre-creonly', dest='sscco', type=str, action='store', required=True, help='output file from standard cre refinement, cre shuffling, and only cre features used in prediciton')
	parser.add_argument('--standard-shuffletss-tssonly', dest='sstto', type=str, action='store', required=True, help='output file from standard cre refinement, tss shuffling, and only tss features used in prediction')

	parser.add_argument('--nearest-control-both', dest='ncb', type=str, action='store', required=True, help='output file from nearest gene assignment, no shuffling, and both features used in prediction')
	parser.add_argument('--nearest-shufflecre-both', dest='nscb', type=str, action='store', required=True, help='output file from nearest gene assignment, cre shuffling, and both features used in prediction')
	parser.add_argument('--nearest-shuffletss-both', dest='nstb', type=str, action='store', required=True, help='output file from nearest gene assignment, tss shuffling, and both features used in prediction')
	parser.add_argument('--nearest-control-creonly', dest='ncco', type=str, action='store', required=True, help='output file from snearest gene assignment, no shuffling, and only cre features used in prediction')
	parser.add_argument('--nearest-control-tssonly', dest='ncto', type=str, action='store', required=True, help='output file from nearest gene assignment, no shuffling, and only tss features used in prediction')
	parser.add_argument('--nearest-shufflecre-creonly', dest='nscco', type=str, action='store', required=True, help='output file from nearest gene assignment, cre shuffling, and only cre features used in prediciton')
	parser.add_argument('--nearest-shuffletss-tssonly', dest='nstto', type=str, action='store', required=True, help='output file from nearest gene assignment, tss shuffling, and only tss features used in prediction')

	parser.add_argument('--outfile', dest='outfile', type=str, action='store', required=True, help='path to store concatenated output file')
	return parser

def add_columns(filename, filetype):
	dtypes = [("Genome", "<U4"), ("CT", "<U30"), ("R2", np.float32)]
	indf = np.loadtxt(filename, dtype=np.dtype(dtypes), skiprows=1)
	filetypewords = {"scb": {"exp": "standard",
							 "condition": "control",
							 "features": "both"},
                    "sscb": {"exp": "standard"},
							  "condition": "shufflecre",
							  "features": "both"},
					"sstb": {"exp": "standard",
							  "condition": "shuffletss",
							  "features": "both"},
                    "scco": {"exp": "standard",
							  "condition": "control",
							  "features": "cre"},
                    "scto": {"exp": "standard",
							  "condition": "control",
							  "features": "promoter"},
                    "sscco": {"exp": "standard",
							  "condition": "shufflecre",
							  "features": "cre"},
                    "sstto": {"exp": "standard",
							  "condition": "shuffletss",
							  "features": "promoter"},
                    "ncb": {"exp": "nearest",
							 "condition": "control",
							 "features": "both"},
                    "nscb": {"exp": "nearest"},
							  "condition": "shufflecre",
							  "features": "both"},
                    "nstb": {"exp": "nearest",
							  "condition": "shuffletss",
							  "features": "both"},
                    "ncco": {"exp": "nearest",
							  "condition": "control",
							  "features": "cre"},
                    "ncto": {"exp": "nearest",
							  "condition": "control",
							  "features": "promoter"},
                    "nscco": {"exp": "nearest",
							   "condition": "shufflecre",
							   "features": "cre"},
                    "nstto": {"exp": "nearest",
							   "condition": "shuffletss",
							   "features": "promoter"}}

	indf['Experiment'] = filetypewords[filetype]["exp"]
	indf['Condition'] = filetypewords[filetype]["condition"]
	indf['Features'] = filetypewords[filetype]["features"]
	indf['Rep'] = 0
	return indf

def concatenate_together(args):
	dtypes = [("Genome", "<U4"), ("CT", "<U30"), ("R2", np.float32),
			  ('Experiment', "<U8" )('Condition', "<U10"), ('Features', "<U8"), ("Rep", np.int32)]
	newarr = pd.DataFrame(dtypes = np.dtype(dtypes))
	for filename, filetype in zip([args.scb, args.sscb, args.sstb, args.scco, args.scto, args.sscco, args.sstto, args.ncb, args.nscb, args.nstb, args.ncco, args.ncto, args.nscco, args.nstto],
									["scb", "sscb", "sstb", "scco", "scto", "sscco", "sstto", "ncb", "nscb", "nstb", "ncco", "ncto", "nscco", "nstto"]):
		new_arr = pd.concat([new_arr, add_columns(filename, filetype)])
	return new_arr

main()
