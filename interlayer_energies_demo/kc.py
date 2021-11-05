import ase
import ase.calculators.lj
from interlayer_energies_demo import generate_geometry
import scipy.optimize
import numpy as np
import numpy.linalg as la
import os

from ase.calculators.lammpsrun import LAMMPS
import load

def get_params_hybrid():
    parameters = {
        'pair_style': 'hybrid/overlay rebo kolmogorov/crespi/full 16.0 1',
        'pair_coeff': [
            '* * rebo CH.rebo C C',
            '* * kolmogorov/crespi/full CH_taper.KC C C'
            ],
        'atom_style': 'full',
        'specorder': ['C', 'C'],
        }
    files = ['CH.rebo', 'CH_taper.KC']
    return parameters, files

def get_params_kc_only():
    parameters = {
        'pair_style': 'hybrid/overlay kolmogorov/crespi/full 16.0 1',
        'pair_coeff': [
            '* * kolmogorov/crespi/full CH_taper.KC C C'
            ],
        'atom_style': 'full',
        'specorder': ['C', 'C'],
        }
    files = ['CH_taper.KC']
    return parameters, files

def write_kc_potential(z0, C0, C2, C4, C, delta, lamda, A):
    with open('CH_taper.KC', 'w') as f:
        f.write(f'''
# Refined parameters for Kolmogorov-Crespi Potential with taper function
#
# Cite as W. Ouyang, D. Mandelli, M. Urbakh and O. Hod, Nano Letters 18, 6009-6016 (2018).
#
#        z0         C0          C2             C4            C            delta        lambda       A           S    rcut
C C  {z0}   {C0}    {C2}     {C4}   {C}    {delta}    {lamda}   {A}    1.0    2.0
''')

e1 = 0

def eval_energy(df, z0, C0, C2, C4, C, delta, lamda, A):
    """
    Finds the energy for a given geometry (in `df`) and KC paramters
    """
    global e1
    os.system('rm -r tmp/')

    energy = []
    for it, row in df.iterrows():
        atoms = generate_geometry.create_graphene_geom(row['d'], row['disregistry'])
        atoms.set_array('mol-id', np.array([0, 0, 1, 1]))
        parameters, files = get_params_kc_only()
        write_kc_potential(z0, C0, C2, C4, C, delta, lamda, A)
        atoms.calc = LAMMPS(files=files, keep_tmp_files=True, tmp_dir='tmp/', **parameters)
        e = atoms.get_potential_energy()/len(atoms)
        energy.append(e)

    en =  np.asarray(energy) - np.min(energy) + np.min(df['energy'])
    error = la.norm(en - e1)
    print(f'energy diff: {error}')
    e1 = en.copy()
    print('params: ', z0, C0, C2, C4, C, delta, lamda, A)
    print('*'*50)
    return en

def fit_kc(df, p0):
    ydata = df['energy']
    popt, pcov = scipy.optimize.curve_fit(eval_energy, df, ydata, p0=p0, method='trf')
    print(popt)
    return popt

# z0_0, C0_0, C2_0, C4_0, C_0, delta_0, lamda_0, A_0
# df = load.load_data()
# eval_energy(df, 3.416084, 20.021583, 10.9055107, 4.2756354, 1.0010836E-2, 0.8447122, 2.9360584, 14.3132588)

def main():
    import load
    df = load.load_data()
    p0 = (3.416084, 20.021583, 10.9055107, 4.2756354, 1.0010836E-2, 0.8447122, 2.9360584, 14.3132588)
    popt = fit_kc(df, p0)
    print(popt)
    return
    # df['en'] = eval_energy(df, *popt)
    # import matplotlib.pyplot as plt
    # import seaborn as sns
    # g = sns.FacetGrid(hue='disregistry', col='disregistry',data =df)
    # g.map( plt.plot,'d','en')
    # g.map( plt.errorbar,'d','energy', 'energy_err', marker='o', mew=1, mec='k', linestyle="")
    # print(df)
    # g.add_legend()
    # plt.show()

if __name__ == '__main__':
    main()
