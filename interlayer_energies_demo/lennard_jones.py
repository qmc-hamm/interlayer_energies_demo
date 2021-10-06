import ase
import ase.calculators.lj
from interlayer_energies_demo import generate_geometry
import scipy.optimize
import numpy as np

def eval_energy(df, sigma, epsilon):
    """ 
    """
    print(sigma, epsilon)
    energy = []
    for it, row in df.iterrows():
        atoms = generate_geometry.create_graphene_geom(row['d'], row['disregistry'])
        calc = ase.calculators.lj.LennardJones(sigma=sigma, epsilon=epsilon)
        atoms.calc=calc
        energy.append(atoms.get_potential_energy()/len(atoms))

    lj_en =  np.asarray(energy)- np.min(energy) + np.min(df['energy'])
    return lj_en
        


def fit_lj(df, sigma0, epsilon0):
    ydata = df['energy']
    popt, pcov = scipy.optimize.curve_fit(eval_energy, df, ydata, p0=(sigma0, epsilon0))
    return popt
    
    


if __name__=="__main__":
    import load
    df = load.load_data()
    sigma, epsilon = fit_lj(df, sigma0=3.5, epsilon0=3e-2)

    df['lj_en'] = eval_energy(df, sigma=sigma, epsilon=epsilon)
    import matplotlib.pyplot as plt
    import seaborn as sns
    g = sns.FacetGrid(hue='disregistry', col='disregistry',data =df)
    g.map( plt.plot,'d','lj_en')
    g.map( plt.errorbar,'d','energy', 'energy_err', marker='o', mew=1, mec='k', linestyle="")
    print(df)
    g.add_legend()
    plt.show()