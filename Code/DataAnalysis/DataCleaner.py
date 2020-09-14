import pandas as pd

# static variables
tool1 = 'Polyphen2'
tool2 = 'PROVEAN'
tool3 = 'SIFT'


# clean all the redundant whitespace in the generated file, then output it
def clean_pph2_data():
    with open('pph2-full.txt', 'r') as file:
        for line in file:
            if line.startswith('##'):  # ignore the comments
                continue
            res = []
            values = line.strip().split('\t')
            for value in values:
                res.append(value.strip())
            f_name = open('pph2-cleaned.txt', 'a+')  # use 'a+' to avoid file overwrite
            print('\t'.join(x for x in res), file=f_name)
    print('done clean pph2')


# load data files and sort them out
def sort_data():
    # load Polyphen2 data
    pph2_col_names = ['#o_acc', 'o_pos', 'o_aa1', 'o_aa2', 'pph2_class']
    pph2 = pd.read_table('pph2-cleaned.txt', sep='\t', usecols=pph2_col_names)
    pph2.columns = ['protein', 'position', 'origin', 'mutated', tool1]
    pph2[tool1] = pph2[tool1].str.capitalize()

    # load SIFT and PROVEAN data
    col_names = ['PROTEIN_ID', 'POSITION', 'RESIDUE_REF',
                 'RESIDUE_ALT', 'PREDICTION (cutoff=-2.5)', 'PREDICTION (cutoff=0.05)']
    # the PROVEAN batch process only allow files smaller than 1MB to be uploaded,
    # so the data set was divided into four
    table1 = pd.read_csv("DataAnalysis/provean1.tsv", sep='\t', usecols=col_names)
    table2 = pd.read_csv("DataAnalysis/provean2.tsv", sep='\t', usecols=col_names)
    table3 = pd.read_csv("DataAnalysis/provean3.tsv", sep='\t', usecols=col_names)
    table4 = pd.read_csv("DataAnalysis/provean4.tsv", sep='\t', usecols=col_names)
    sift_and_provean = pd.concat([table1, table2, table3, table4])
    sift_and_provean.columns = ['protein', 'position', 'origin', 'mutated', tool2, tool3]

    # Data combination
    data_set = pd.merge(pph2, sift_and_provean, on=['protein', 'position', 'origin', 'mutated'])

    # for SIFT, 'damaging' is the same as 'deleterious' and 'tolerated' is the same as 'neutral'
    data_set.loc[data_set[tool3] == 'Damaging', tool3] = 'Deleterious'
    data_set.loc[data_set[tool3] == 'Tolerated', tool3] = 'Neutral'

    # add z-score to the data
    protein_z = {'P51587': 'random', 'P40692': 'random', 'P38398': 'random', 'P43246': 'random', 'P52701': 'random',
                 'Q6P2Q9': 'high', 'Q07954': 'high', 'Q9Y4A5': 'high', 'Q00610': 'high',
                 'O75841': 'low', 'Q3KR16': 'low', 'Q96LZ3': 'low', 'Q9Y6M9': 'low', 'O43823': 'low'}
    data_set['zscore'] = data_set['protein'].apply(lambda x: protein_z.get(x))

    # export the data to a file for further use
    data_set.to_csv("dataSet.csv", index=False)
    print('export success')


if __name__ == '__main__':
    clean_pph2_data()
    sort_data()
