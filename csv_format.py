import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import ticker

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
    parameters = {'upstream': range(-1000000, 1000000, 1000), 'downstream': range(-1000000, 1000000, 1000), 'height': [0,0.25, 0.5,0.75, 1.0,0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0 ,4.5, 5.0] + range(5, 301)}
    # parameters = {'upstream': range(-1000000, 1000000, 1000), 'downstream': range(-1000000, 1000000, 1000),
    #               'height': range(5, 301)}
    markers = ['h3k4me1', 'h3k27ac', 'h3k4me3', 'h3k27me3']
    # markers = ['h3k4me1', 'h3k4me3', 'h3k27ac']
    # markers = ['h3k27me3']
    features = ['total_width', 'height', 'total_signal', 'kurtosis', 'skewness']
#['total_width', 'single_width',  'total_signal', 'kurtosis', 'skewness', 'single_signal',]
    for marker in markers:
        for feature in features:
            for parameter in parameters.keys():
                results = []
                df = pd.read_csv(directory + 'grid_path_'+marker+'_'+feature+'.csv')
                for p in parameters[parameter]:
                    cur_df = df[df[parameter] == p]
                    for other_parameter in parameters.keys():
                        if other_parameter != parameter:
                            # print parameters[other_parameter]
                            cur_df = cur_df[cur_df[other_parameter].isin(parameters[other_parameter])]
                    if cur_df.shape[0] == 0:
                        continue
                    results.append((p, cur_df['logP'].min()*-1))
                df = pd.DataFrame(results)
                df.columns = [parameter, 'best_logP']
                df.to_csv(directory+marker+'_'+feature+'_'+parameter+'.csv', index=None)

def get_best_parameter(directory='./csv/'):
    results = {}
    paths = [x for x in os.listdir(directory) if x.find('change') == -1 and x.endswith('.csv')]
    for path in paths:
        df = pd.read_csv(directory+path)
        marker = path.split('_')[2]
        info = path[:-4].split('_')
        feature = '_'.join(info[3:])
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

def generate_plot(df_path, marker, verbose=True):
    df = pd.read_excel(df_path)
    df = df[df[marker]!=0]
    if marker == 'upstream':
        df[marker] = df[marker] * -1
    df = df.set_index([marker])
    columns = df.columns
    colors = ['red', 'blue', 'green', 'purple', 'yellow']
    shapes = ['o', 'D','^', 'x']
    labels = []
    handals = []

    fig, ax = plt.subplots()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    if marker != 'height':
        ax.set_xscale('symlog')
        ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))

    for i in range(0, len(columns)):
        shape = shapes[i]
        # ax = df.plot.scatter(columns[0], columns[i], color=colors[i], xlim=(df[columns[0]].min(), df[columns[0]].max()),
        #                      ylim=(0, df[columns[i]].max()),
        #                       edgecolor='', ax=ax)
        if marker == 'downstream':
            ax = plt.scatter(df.index/1000, df[columns[i]], color=colors[i], edgecolor='', marker=shape)
            mask = np.isfinite(df[columns[i]])
            plt.plot(df.index[mask]/1000, df[columns[i]][mask], color=colors[i], linestyle='-')

        elif marker == 'upstream':
            ax = plt.scatter(df.index / 1000, df[columns[i]], color=colors[i], edgecolor='', marker=shape)
            mask = np.isfinite(df[columns[i]])
            plt.plot(df.index[mask]/1000, df[columns[i]][mask], color=colors[i], linestyle='-')
        else:
            ax = plt.scatter(df.index, df[columns[i]], color=colors[i], edgecolor='', marker=shape)
            mask = np.isfinite(df[columns[i]])
            plt.plot(df.index[mask], df[columns[i]][mask], color=colors[i], linestyle='-')

        handals.append(ax)
        labels.append(columns[i])

    # print handles, labels
    # ax.legend(handles, labels, loc='best')
    font = {'fontname': 'Helvetica','fontsize':10}

    if marker == 'downstream':
        plt.xlabel('Downstream distance (kb)', **font)
        plt.xlim((df.index.min()/1000, df.index.max()/1000))
    elif marker == 'upstream':
        plt.xlabel('Upstream distance (kb)', **font)
        plt.xlim(df.index.max() / 1000, df.index.min() / 1000)
    else:
        plt.xlabel('Height cutoff', **font)
        # plt.xlim((np.log10(df.index.min()), np.log10(df.index.max())))
        # plt.xlim((0, df.index.max()))
        plt.xlim((0, 100))
    #
    plt.ylabel('-log10 enrich P', **font)
    # print labels
    # print ((0, df.max().max()))
    plt.ylim(((0, (int(df.max().max())/10+1)*10)))
    plt.xticks(**font)
    plt.yticks(**font)
    print labels
    if verbose:
        plt.legend(handals, labels, scatterpoints=1, loc=9, fontsize=10, ncol=2)

    plt.savefig(df_path.replace('.xlsx', '.pdf'))
    plt.close('all')

def group_plot(marker, ranges, path):
    csvs = [x for x in os.listdir(path) if x.find(marker)!= -1 and not x.startswith('grid') and not x.endswith('.xlsx') and x.endswith('.csv')]
    # print csvs
    csvs = [x for x in csvs if x.split('_')[2].find(marker)!=-1]
    # print csvs, marker
    dfs = [pd.read_csv(path+'/'+ c) for c in csvs]
    # names = [c.split('_')[0]+'_'+''.join(c.split('_')[2:]).replace('onco', 'oncogene').replace('sup', 'suppressor').replace('inputgenebody', '_genebody').replace('input', "_TSS").replace(marker, '').replace('.csv','').replace('width','') for c in csvs]
    names = [c.split('_')[0] for c in csvs]
    print names
    for i in range(len(dfs)):
        dfs[i] = dfs[i].drop_duplicates()
        dfs[i] = dfs[i].set_index([marker])
    results = []
    for r in ranges:
        cur_results = [r]
        for i in range(len(dfs)):
            try:
                cur_results.append(dfs[i].ix[r, 'best_logP'])
            except:
                cur_results.append(None)
        results.append(cur_results)
    result_df = pd.DataFrame(results)
    result_df.columns=[marker] + names
    result_df.to_excel(path+'/'+marker+'.xlsx', index=False)

    generate_plot(path+'/'+marker+'.xlsx', marker)
    return





# reformat()
#
#
# get_best()
# get_best_parameter()
#
#
# group_plot('height', [0,0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9,0.25,0.5,0.75]+ range(1,50), './csv/h3k4me3_input_oncogene_tumor/')
# group_plot('height', range(1,50), './csv/mid/height/')
# group_plot('downstream', range(0,1000000,1000), './csv/mid/height')
# group_plot('upstream', range(-1000000,0, 1000), './csv/mid/height')

# group_plot('height', [0,0.1,0.2,0.3,0.4,0.6,0.7,0.8,0.9,0.25,0.5,0.75]+ range(1,200), './csv/genebody/single_signal/')
group_plot('height', range(1,100), './csv/6th/kurtosis/')
group_plot('downstream', range(0,500000,1000), './csv/6th/kurtosis/')
group_plot('upstream', range(-500000,0, 1000), './csv/6th/kurtosis/')
