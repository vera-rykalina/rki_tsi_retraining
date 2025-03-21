#!/usr/bin/env python3

# Import libraries
import argparse
import pandas as pd


# Global variables 
_args = None

def initialise():
    '''
    Parse command-line arguments.
    '''
    global _args
    parser = argparse.ArgumentParser( description="Concatinate MAF files for individual samples" )
    parser.add_argument( "-r", "--reference", required=True, help="A comma-separated XHB2 reference file" )
    parser.add_argument( "-l", "--mafs", nargs='+', required=True, help="A list of comma-separated calculated MAF files" )
    parser.add_argument( "-o", "--outpath", required=True, help="Filename for saving joined MAFs." )
    _args = parser.parse_args()
    return

def load_xhb2_reference(fpath):
    ''' Load reference. '''
    ref_df = pd.read_csv(fpath, sep=",")
    print("Loaded HXB2 reference, shape={}".format(ref_df.shape))
    return ref_df


def loads_mafs(fpath):
    ''' Load mafs. '''
    # Create an empty list of dfs with mafs
    maf_dfs = []
    for maf in fpath[0:]:
        maf_dfs.append(pd.read_csv(maf, sep = ",", index_col=False))
        print(maf.head())
        print(maf.shape)
    print(maf_dfs)
    return maf_dfs