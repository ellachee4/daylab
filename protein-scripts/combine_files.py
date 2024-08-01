import pandas as pd

df_antibodies = pd.read_csv('antibody_availability_data.csv')
df_go = pd.read_csv('protein_article_counts2.csv')

df_all = pd.concat([df_antibodies, df_go], axis = 1)
df_all.drop(columns=['UNIPROT'], inplace=True)

for col in ['Number of Antibodies', 'Article Count', 'Interactions', 'GO Score']:
    df_all[f'{col} (normalized)'] = df_all[col]/df_all[col].std()

def compute_score(row):
    score = row['Article Count (normalized)']*-2 + row['Number of Antibodies (normalized)']*0.5 + row['Interactions (normalized)'] + row['GO Score (normalized)']
    return score

df_all['Overall Score'] = df_all.apply(lambda x: compute_score(x), axis=1)
cols = ['Uniprot', 'Protein Symbol','Antibody Link', 'Number of Antibodies', 'Number of Providers',
        'Article Count', 'Interactions', 'GO Score', 'Overall Score', 'GO Terms']

df_all = df_all[cols]
df_all.to_csv('combined-scores.csv')