import pandas as pd
import statsmodels.api as sm
import numpy as np

# static variables
tool1 = 'Polyphen2'
tool2 = 'PROVEAN'
tool3 = 'SIFT'

# Dictionary for ORF detection
molecular_weight = {'A': 71.08, 'R': 156.19, 'N': 114.11, 'D': 115.09, 'C': 103.15, 'E': 129.12,
                    'Q': 128.13, 'G': 57.05, 'H': 137.14, 'O': 113.11, 'I': 113.16, 'L': 113.16,
                    'K': 128.18, 'M': 131.20, 'F': 147.18, 'P': 97.12, 'U': 121.09, 'S': 87.08,
                    'T': 101.11, 'W': 186.22, 'Y': 163.18, 'V': 99.13}

# protein z score
protein_z_score = {'P51587': -1.29, 'P40692': -0.304, 'P38398': 0.582, 'P43246': -2.45, 'P52701': -2.78,
                   'Q6P2Q9': 8.28, 'Q07954': 8.25, 'Q9Y4A5': 8.17, 'Q00610': 7.76,
                   'O75841': 0.000674, 'Q3KR16': 0.000607, 'Q96LZ3': 0.000593, 'Q9Y6M9': 0.000386, 'O43823': 0.0000295}


def create_data_set(tool_name):
    # import data
    df = pd.read_csv("DataAnalysis/dataSet.csv", usecols=[tool_name, 'position', 'protein', 'origin', 'mutated'])
    df[tool_name] = df.apply(lambda x: int(x[tool_name] != 'Neutral'), axis=1)
    df['origin'] = df.apply(lambda x: molecular_weight[x.origin], axis=1)
    df['mutated'] = df.apply(lambda x: molecular_weight[x.mutated], axis=1)
    df['z_score'] = df.apply(lambda x: protein_z_score[x.protein], axis=1)
    df['intercept'] = 1.0
    return df[[tool_name, 'origin', 'position', 'mutated', 'z_score', 'intercept']]


def get_trained_model(df, tool_name):
    train_cols = df.columns[1:]
    logistic = sm.Logit(df[tool_name], df[train_cols])
    return logistic.fit()


def predict(model, df):
    predict_cols = df.columns[1:]
    df['predict'] = model.predict(df[predict_cols])
    return df


def evaluate_prediction(df):
    size = 0
    n = 0
    for value in df.values:
        res = value[-1]  # the result created by the model, stored in the last column
        real = int(value[0])  # the real result
        # if the probability is over the function value then the prediction succeeded
        size += 1
        if res > 0.66 and real == 1:
            n += 1
        elif res < 0.66 and real == 0:
            n += 1
    return size, n


if __name__ == "__main__":
    # get the selected tool's data
    data = create_data_set(tool3)

    tdf = data[0:-1:5]  # training data
    cdf = data[0:-1:3]  # test data
    print(data.head())

    # train the model
    model = get_trained_model(tdf, tool3)

    # make prediction
    cdf = predict(model, cdf)

    # get the evaluation
    total, hit = evaluate_prediction(cdf)

    # output
    print('Total: %d, Hit: %d, Precision: %.2f' % (total, hit, 100.0 * hit / total))

    # Total: 49, Hit: 30, Precision: 61.22
    print(model.summary())

    params = model.params
    conf = model.conf_int()
    conf['OR'] = params
    conf.columns = ['2.5%', '97.5%', 'OR']
    print(np.exp(conf))
