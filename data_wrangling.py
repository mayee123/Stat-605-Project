import pandas as pd

#read in csv files
df1 = pd.read_csv('typing.csv')
df2 = pd.read_csv('nuccore.csv')
df3=pd.read_csv('amr.tsv',sep='\t')

#See number of each observed host range taxonomic rank
unique_values = df1['observed_host_range_ncbi_rank'].value_counts()
print(unique_values)

#create new dataframe called df_amr with the number of Antimicrobial Resistance Genes
amr=df3['NUCCORE_ACC'].value_counts()
df_amr=(pd.DataFrame({'NUCCORE_ACC':amr.index, 'arg_count':amr.values}))

#merge the dataframes on the 'NUCCORE_ACC' row, filling in 0 for ARG_COUNT if the accession does not have any Antimicrobial Resistance Genes
merged_df1 = pd.merge(df1, df2, on='NUCCORE_ACC')
merged_df2 = pd.merge(merged_df1,df_amr, on='NUCCORE_ACC',how='left')
merged_df2['arg_count'] = merged_df2['arg_count'].fillna(0)

#filter out any accessions without 'observed_host_range_ncbi_rank'
df_filtered = merged_df2[merged_df2['observed_host_range_ncbi_rank'].notnull()]

#create new column 'range' with value 0 if observed_host_range_ncbi_rank is 'genus' and 1 otherwise
df_filtered['range'] = df_filtered['observed_host_range_ncbi_rank'].apply(
    lambda x: 0 if x == 'genus' else 1)

#create csv file called 'plasmid_data.csv'
df_filtered.to_csv('plasmid_data.csv', index=False)

