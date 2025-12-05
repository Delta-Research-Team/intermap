# Created by rglez at 7/17/25
"""
This script is intended to process the InterMap output pickle of the channel5
simulation to report average prevalence of interactions either in compact
or extended format.
"""
import pandas as pd
from bitarray import bitarray as ba

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
    df_temp = df.iloc[value]
    df_temp['s1'] = df_temp['s1'].str.split('_').str[:2].str.join('_')
    df_temp['s2'] = df_temp['s2'].str.split('_').str[:2].str.join('_')
    df_temp = df_temp[~df_temp['inter_name'].isin(to_exclude)]

    compacted_info = []
    groups_s1_s2 = df_temp.groupby(['s1', 's2']).indices
    for (s1, s2), indices in groups_s1_s2.items():
        timeseries = df_temp.iloc[indices]['time'].apply(ba).tolist()
        joint = timeseries[0]
        for bit_time in timeseries[1:]:
            joint |= bit_time
        prevalence = joint.count() / len(joint) * 100
        compacted_info.append((s1, s2, prevalence))
    joint_df = pd.DataFrame(compacted_info,
                            columns=['s1', 's2', f'prev-{key}'])
    dfs.append(joint_df)

    for _, row in joint_df.iterrows():
        all_combinations.add((row['s1'], row['s2']))

# Merge all dataframes and calculate the average prevalence
merged_df = pd.DataFrame(list(all_combinations),
                         columns=['s1', 's2'])

for df_group in dfs:
    merged_df = merged_df.merge(
        df_group,
        on=['s1', 's2'],
        how='left'
    )
merged_df.fillna(0, inplace=True)
merged_df['prev-mean'] = merged_df.filter(like='prev-').mean(axis=1)

# Sort and exclude some values of the inter_name column
merged_df['resid'] = merged_df['s1'].str.split('_', expand=True)[1].apply(int)

merged_df.sort_values(by=['resid', 'prev-mean'], inplace=True,
                      ascending=[True, False])

# Save the results to a CSV file
output_file = f'{output_dir}/channel5_compact.csv'
merged_df.to_csv(output_file, index=False)
