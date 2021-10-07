import ase
import ase.calculators.lj
from interlayer_energies_demo import generate_geometry
import scipy.optimize
import numpy as np
import os

from ase.calculators.lammpsrun import LAMMPS
import load

def eval_energy(df, z0, C0, C2, C4, C, delta, lamda, A):
    """
    """

    os.system('rm -r tmp_dir/')
    energy = []
    for it, row in df.iterrows():
        atoms = generate_geometry.create_graphene_geom(row['d'], row['disregistry'])
        print('\n'.join(dir(atoms)))
        parameters = {
            'pair_style': 'hybrid/overlay rebo kolmogorov/crespi/full 16.0 1',
            'pair_coeff': [
                # '* * none',
                '* * rebo CH.rebo C C',
                '* * kolmogorov/crespi/full CH_taper.KC C C'
                ],
            'atom_style': 'full',
            'specorder': ['C', 'C'],
            'mol-id': [1, 1, 2, 2]
            }
        files = ['CH.rebo', 'CH_taper.KC']
        lammps = LAMMPS(parameters=parameters, files=files, keep_tmp_files=True, tmp_dir='tmp_dir/')
        atoms.calc = lammps
        print('*'*50)
        print(atoms.get_potential_energy()/len(atoms))
        print('after')
        energy.append(atoms.get_potential_energy()/len(atoms))
        break

    lj_en =  np.asarray(energy) - np.min(energy) + np.min(df['energy'])
    return lj_en


df = load.load_data()
eval_energy(df, 3.416084, 20.021583, 10.9055107, 4.2756354, 1.0010836E-2, 0.8447122, 2.9360584, 14.3132588)
