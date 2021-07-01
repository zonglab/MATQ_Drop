#!/usr/local/bin/python3.5

import sys
import pandas as pd

UMI_matrix, BC_read_count, Output_File_Name = sys.argv[1], sys.argv[2], sys.argv[3]

Read_count=pd.read_csv(BC_read_count,sep="\t", index_col=0, header=None).T.squeeze()

df = pd.read_csv(UMI_matrix, sep="\t", index_col=0)

UMI_count = df.fillna(0).sum()
Gene_count = df.fillna(0).astype(bool).sum()

a = pd.concat([Read_count.astype(int),UMI_count.astype(int),Gene_count.astype(int)],axis=1)
a = a.fillna(0).astype(int)
a.columns = ['Raw_Reads', 'UMI_count', 'Gene_count']
a.to_csv(Output_File_Name, sep='\t')


