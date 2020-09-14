import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# static variables
tool1 = 'Polyphen2'
tool2 = 'PROVEAN'
tool3 = 'SIFT'


# perform an overall analysis on the data
def plot_aly_overall(data):
    del data['protein']
    # reformat the data matrix
    pp2 = data.groupby('zscore')[tool1].value_counts().unstack().reset_index()
    pp2['tool'] = tool1
    prv = data.groupby('zscore')[tool2].value_counts().unstack().reset_index()
    prv['tool'] = tool2
    sft = data.groupby('zscore')[tool3].value_counts().unstack().reset_index()
    sft['tool'] = tool3

    # merge three data frames
    data = pd.concat([pp2, prv, sft])
    data.reset_index(inplace=True)
    del data['index']
    print(data)
    data['Ratio'] = data.apply(lambda x: x.Deleterious / (x.Deleterious + x.Neutral), axis=1)
    del data['Deleterious']
    del data['Neutral']
    data = data.groupby(['zscore', 'tool']).sum().reset_index()
    data = data.groupby(['zscore', 'tool'], sort=True).sum()
    data = data.unstack(1)
    data.columns = [tool1, tool2, tool3]
    data.columns.names = ['Tool']
    print(data)
    data.plot(kind='bar', figsize=(10, 6), title='The Ratio of Pathogenic Prediction', stacked=False)
    plt.show()


# the scoring system evaluates each tool based on their results comparing to other two tools'
def score_system(a, b, c):
    # the format of the scoring log is as follows:
    #          Tool_1  Tool_2  Tool_3
    #  Match       t1      t2      t3
    #  Negative    s1      s2      s3
    #  Positive    l1      l2      l3
    scoring_log = np.zeros(9)
    if a == b == c:  # we believe the result is trustworthy
        scoring_log[0] += 1
        scoring_log[1] += 1
        scoring_log[2] += 1
    else:
        if a == 'Deleterious':
            scoring_log[3] += 1
        else:
            scoring_log[6] += 1
        if b == 'Deleterious':
            scoring_log[4] += 1
        else:
            scoring_log[7] += 1
        if c == 'Deleterious':
            scoring_log[5] += 1
        else:
            scoring_log[8] += 1
    return scoring_log


# apply the scoring system to the tools
def scoring(data):
    x = data[tool1].values
    y = data[tool2].values
    z = data[tool3].values
    res = np.zeros(9)
    for i in range(len(x)):
        res += score_system(x[i], y[i], z[i])
    res = pd.DataFrame(np.array(res).reshape(3, 3), columns=[tool1, tool2, tool3])
    res.index = ['Match', 'Neutral', 'Deleterious']
    res['t1'] = res.apply(lambda dat: dat.Polyphen2 / sum(res[tool1]), axis=1)
    res['t2'] = res.apply(lambda dat: dat.PROVEAN / sum(res[tool2]), axis=1)
    res['t3'] = res.apply(lambda dat: dat.SIFT / sum(res[tool3]), axis=1)
    del res[tool1]
    del res[tool2]
    del res[tool3]
    return res


def plot_scoring(data):
    # plot the overview of data
    data = data.groupby(['correctness', 'zscore'], sort=True).sum()[[tool1, tool2, tool3]].reset_index(
        level=1, drop=False)
    zs = data.zscore.unique().tolist()
    fig, axes = plt.subplots(ncols=len(zs), nrows=1)
    for z, ax in zip(zs, axes.ravel()):
        data.loc[data['zscore'] == z].plot(ax=ax, kind='bar', sharey=True, ylim=[0, 0.8], figsize=(10, 6), title=z)
    plt.tight_layout()
    plt.subplots_adjust(top=0.88)
    plt.suptitle('The ratio of each tool having non-pathogenic and pathogenic prediction',)
    plt.show()


if __name__ == "__main__":

    # perform an overall analysis on the data
    df = pd.read_csv("dataSet.csv", usecols=['protein', tool1, tool2, tool3, 'zscore'])
    plot_aly_overall(df)

    # perform the scoring system on the data
    scoring_data = df.groupby('zscore').apply(lambda data: scoring(data)).reset_index()
    scoring_data.columns = ['zscore', 'correctness', tool1, tool2, tool3]
    plot_scoring(scoring_data)
    print(scoring_data)
