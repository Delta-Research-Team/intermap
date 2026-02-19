# Created by rglez at 7/17/25
"""
This script is intended to process the InterMap output pickle of the channel5
simulation to report average prevalence of interactions either in compact
or extended format.
"""
import pandas as pd

# =============================================================================
# User defined variables
# =============================================================================
csv = '/media/rglez/Roy2TB/Dropbox/RoyData/intermap/tutorial-mayank/outputs/ligs-channel_InterMap.csv'
output_dir = '/media/rglez/Roy2TB/Dropbox/RoyData/intermap/tutorial-mayank/outputs/'
to_exclude = []
# =============================================================================

# Separate the csv by note1 (segname of the ligand)
df = pd.read_csv(csv)
grouped = df.groupby(['note1']).indices
dfs = []
all_combinations = set()
for key, value in grouped.items():
    df_temp = df.iloc[value][['note1', 's1', 's2', 'inter_name', 'prevalence']]
    df_temp['s1'] = df_temp['s1'].str.split('_').str[:2].str.join('_')
    df_temp['s2'] = df_temp['s2'].str.split('_').str[:2].str.join('_')
    df_id = df_temp.iloc[0, 0]
    df_temp.rename(columns={'prevalence': f'prev-{df_id}'}, inplace=True)
    dfs.append(df_temp[['s1', 's2', 'inter_name', f'prev-{df_id}']])

    for _, row in df_temp.iterrows():
        all_combinations.add((row['s1'], row['s2'], row['inter_name']))

# Merge all dataframes and calculate the average prevalence
merged_df = pd.DataFrame(list(all_combinations),
                         columns=['s1', 's2', 'inter_name'])
for df_group in dfs:
    merged_df = merged_df.merge(
        df_group,
        on=['s1', 's2', 'inter_name'],
        how='left'
    )
merged_df['prev-mean'] = merged_df.filter(like='prev-').mean(axis=1)

# Sort and exclude some values of the inter_name column
merged_df['resid'] = merged_df['s1'].str.split('_', expand=True)[1].apply(int)

merged_df.sort_values(by=['resid', 'prev-mean', 'inter_name'], inplace=True,
                      ascending=[True, False, True])

merged_df = merged_df[~merged_df['inter_name'].isin(to_exclude)]
merged_df.fillna(0, inplace=True)

# Save the results to a CSV file
output_file = f'{output_dir}/channel5_extended.csv'
merged_df.to_csv(output_file, index=False)
