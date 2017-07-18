import os


import pandas as pd
import matplotlib.pyplot as plt

def reformat(directory='./csv/'):
    paths = [x for x in os.listdir(directory) if x.find('change') == -1 and x.endswith('.csv')]
    for path in paths:
        results = []
        print path
        df = pd.read_csv(directory+ path, index_col=0)

        print df

        df.columns = ['p','logP']

        for i in range(df.shape[0]):
           info = df.ix[i, 'p']
           info = info.replace('[', '')
           info = info.replace(']', '')
           cur_result = info.split() + [df.ix[i, 'logP']]
           results.append(cur_result)

        df = pd.DataFrame(results)
        df.columns = ['upstream', 'downstream', 'height', 'logP']
        # df = df[df['logP'] > -82]
        df.to_csv(directory+path, index=None)

def get_best(directory='./csv/'):
    parameters = {'upstream': range(-250000, 250000, 1000), 'downstream': range(-250000, 250000, 1000), 'height': [0.25,0.5,0.75,1.25,1.5,1.75,2.25,2.5,2.75,3.25,3.5,3.75,4.5,4.25,4.75]+range(1, 201, 1)}
    # markers = ['h3k4me1', 'h3k4me3', 'h3k27me3', 'h3k27ac']
    markers = ['h3k27ac']
    features = ['single_width', 'total_signal', 'total_width', 'single_signal']
#['total_width', 'single_width',  'total_signal', 'height', 'kurtosis', 'skewness', 'single_signal',]
    for marker in markers:
        for feature in features:
            for parameter in parameters.keys():
                results = []
                df = pd.read_csv(directory + 'grid_path_'+marker+'_'+feature+'.csv')
                for p in parameters[parameter]:
                    cur_df = df[df[parameter] == p]
                    if cur_df.shape[0] == 0:
                        continue
                    results.append((p, cur_df['logP'].min()*-1))
                df = pd.DataFrame(results)
                df.columns = [parameter, 'best_logP']
                df.to_csv(directory+marker+'_'+feature+'_'+parameter+'_grid_change.csv', index=None)

def get_best_parameter(directory='./csv/'):
    results = {}
    paths = [x for x in os.listdir(directory) if x.find('change') == -1 and x.endswith('.csv')]
    for path in paths:
        df = pd.read_csv(directory+path)
        marker = path.split('_')[2]
        info = path[:-4].split('_')
        if info[-2] == 'single' or info[-2] == 'total':
            feature = info[-2] + '_' + info[-1]
        else:
            feature = info[-1]
        best_p = df['logP'].min()
        cur_df = df[df['logP'] == best_p]
        if (marker, feature) in results:
            if best_p < results[(marker, feature)][-1]:
                results[(marker, feature)] = [tuple(x) for x in cur_df.to_records(index=False)][0]
        else:
            results[(marker, feature)] = [tuple(x) for x in cur_df.to_records(index=False)][0]
    final = []
    for key, value in results.items():
        final.append([key[0], key[1]] + list(results[key]))

    final_df = pd.DataFrame(final)
    final_df.columns = ['marker', 'feature', 'upstream', 'downstream', 'height', 'logP']
    final_df.to_csv('best_parameters_CIG.csv', index=None)

def generate_plot(df_path):
    df = pd.read_csv(df_path)
    columns = df.columns
    df.plot.scatter(columns[0], columns[1], c='black')

    plt.savefig(df_path.replace('.csv', '.png'))

# reformat()
# df1 = pd.read_csv('./csv/1stRound/h3k27ac/grid_path_h3k27ac_single_width.csv')
# df2 = pd.read_csv('./csv/grid_path_h3k27ac_single_width.csv')
#
# df = df1.append(df2)
#
# df.to_csv('./csv/grid_path_h3k27ac_single_width.csv', index=None)


get_best()
get_best_parameter()

# files = ['./csv/1stRound/h3k4me3/' + x for x in os.listdir('./csv/1stRound/h3k4me3/') if x.startswith('h3k4me3') and x.endswith('.csv')]
#
# for f in files:
#     generate_plot(f)