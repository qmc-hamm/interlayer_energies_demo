import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

import load

plt.rc('font', family='serif')

def create_fitted_df():
    df_orig = load.load_data()
    popt = (3.416084, 20.021583, 10.9055107, 4.2756354, 1.0010836E-2, 0.8447122, 2.9360584, 14.3132588)
    df_orig['en'] = eval_energy(df_orig, *popt)
    df_orig['label'] = 'orig'

    # popt = (3.4137495175365897, 19.411451757497794, 11.309475638472934, 6.21195315450154, 0.2863398313122875, 0.8003381594877969, 3.0961781292223756, 14.173639718483972) # results for gtol=1e-04
    df_fit = load.load_data()
    popt = (3.372485976907569, 18.172084448730008, 3.573924286622721, 12.909862639795223, 5.68180602644627, 0.6484292762516102, 3.1279865162981135, 15.244928020134566) # results for gtol=1e-06
    df_fit['en'] = eval_energy(df_fit, *popt)
    df_fit['label'] = 'qmc'
    df = pd.concat([df_orig, df_fit], ignore_index=True)
    df['stacking'] = df['disregistry'].map({0.0: 'AB', 0.16667: 'SP', 0.5: 'Mid', 0.66667: 'AA'})
    print(df)
    df.to_csv('kc.csv', index=False)

def plot():
    df = pd.read_csv('kc.csv')
    df['label'] = df['label'].map({'orig': 'DFT fitted', 'fit': 'QMC fitted'})
    print(df)
    g = sns.FacetGrid(data=df, col='stacking', hue='label', sharey='row')
    g.map(plt.plot, 'd', 'en')
    g.map(plt.errorbar, 'd', 'energy', 'energy_err', marker='o', alpha=0.7, ms=3, linestyle="")
    for i, (stacking, ax) in enumerate(g.axes_dict.items()):
        if i == 0:
            ax.set(ylabel='$E$ (eV/atom)')
    g.set(xlabel='$d$ (ang)')
    g.add_legend(bbox_to_anchor=(0.9, 0.7), fontsize=9)
    plt.savefig('kc_refit.pdf', bbox_inches='tight')

if __name__ == '__main__':
    plot()
