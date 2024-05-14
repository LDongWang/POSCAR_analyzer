class Atom:
    def __init__(self, ele, no, pos):
        self.ele = ele
        self.no = no
        self.pos = pos
        self.name = ele+no
    def __str__(self):
        return f"{self.name}{self.no}:{self.pos}"

    
def adj_dist(Atom1, Atom2, lattice_vector, no_shift=False):
    if no_shift:
        return np.linalg.norm((Atom1.pos-Atom2.pos)@lattice_vector)
    else:
        unit_x = np.array([1.0, 0.0, 0.0])
        x_shift = [-unit_x, 0.0, unit_x]
        unit_y = np.array([0.0, 1.0, 0.0])
        y_shift = [-unit_y, 0.0, unit_y]
        unit_z = np.array([0.0, 0.0, 1.0])
        z_shift = [-unit_z, 0.0, unit_z]
        dist = Atom1.pos-Atom2.pos
        dist_list = []
        for x in x_shift:
            for y in y_shift:
                for z in z_shift:
                    dist_list.append(np.linalg.norm((dist+x+y+z)@lattice_vector))
        return min(dist_list)
    
def lattice_transform(lattice, point_group='p1'):
    '''
    Now support p1 and c2/m symmetry. Other symmetry are under progress.
                                                    2023/09/29 Lidong W.
    '''
    if point_group == 'p1':
        return lattice
    elif point_group == 'c2/m':
        return [lattice[0],2*lattice[1]-lattice[0],lattice[2]]
    else:
        print('Not implemented symmetry - Treat as p1')
        return lattice

def lattice_info(lattice_vector):
    lattice_norm = [np.linalg.norm(i) for i in lattice_vector]
    a,b,c = lattice_vector
    a_norm,b_norm,c_norm = lattice_norm
    alpha = np.rad2deg(np.arccos(b@c/b_norm/c_norm))
    beta = np.rad2deg(np.arccos(a@c/a_norm/c_norm))
    gamma = np.rad2deg(np.arccos(a@b/a_norm/b_norm))
    volume = np.cross(a,b)@c.T
    return [lattice_norm, [alpha, beta, gamma], volume]
    
    
import numpy as np
import argparse
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser('''python POSCAR2distance.py''')
parser.add_argument('filename', help='input POSCAR name')
parser.add_argument('--output', default='POSCARinfo', help='output filename')
parser.add_argument('--symmetry', default='p1', help='specify symmetry point group')
parser.add_argument('--dist', action='store_true', help='output distance matrice')
parser.add_argument('--plot', action='store_true', help='plot bond distribution')
parser.add_argument('--pair', nargs='+', help='Atoms pair to plot, e.g., Ti-O Zr-O')
parser.add_argument('--bandwidth', type=float, default=0.20, help='Bandwidth for KDE. Default is 0.20.')
parser.add_argument('--info', action='store_true', help='output lattice information')
parser.add_argument('--noshift', action='store_true', help='specify no periodic boundary condition')
args = parser.parse_args()
with open(args.filename,'r') as f:
    title = f.readline()
    scale = float(f.readline())
    a = np.array(f.readline().split(),dtype=float)
    b = np.array(f.readline().split(),dtype=float)
    c = np.array(f.readline().split(),dtype=float)
    origin_lattice_vector = np.array([a,b,c])*scale
    lattice_vector = lattice_transform(origin_lattice_vector, args.symmetry)
    ele_names = f.readline().split()
    ele_nums  = np.array(f.readline().split(),dtype=int)
    mode = f.readline().split()[0]
    atm_list = []
    for i in range(len(ele_names)):
        for j in range(ele_nums[i]):
            pos = np.array(f.readline().split(),dtype=float)
            atm_list.append(Atom(ele_names[i],str(j+1),pos))
with open(args.output, mode='w', encoding='utf-8', newline='\n') as f:
    f.write(title+'\n')
    if args.info:
        lattice_norm, lattice_angles, lattice_volume = lattice_info(lattice_vector)
        f.write('Lattice Information\n')
        f.write('cell_length_a:\t{:.3f}\n'.format(lattice_norm[0]))
        f.write('cell_length_b:\t{:.3f}\n'.format(lattice_norm[1]))
        f.write('cell_length_c:\t{:.3f}\n'.format(lattice_norm[2]))
        f.write('cell_angle_alpha:\t{:.3f}\n'.format(lattice_angles[0]))
        f.write('cell_angle_beta:\t{:.3f}\n'.format(lattice_angles[1]))
        f.write('cell_angle_gamma:\t{:.3f}\n'.format(lattice_angles[2]))
        f.write('Lattice_volume:\t{:.3f}\n'.format(lattice_volume))
    if args.dist:
        f.write('Distance Matrix\n')
        f.write(' '.join(ele_names)+'\n')
        f.write(' '.join(str(i) for i in ele_nums)+'\n')
        f.write('\t'+'\t'.join(atm.name for atm in atm_list)+'\n')
        dist_matrix = np.zeros((len(atm_list),len(atm_list)))
        if args.plot:
            dist_list = {}
        for i in range(len(atm_list)):
            left_name = atm_list[i].name
            f.write(left_name+'\t')
            for j in range(len(atm_list)):
                right_name = atm_list[j].name
                dist_matrix[i][j] = adj_dist(atm_list[i],atm_list[j],origin_lattice_vector,no_shift=args.noshift)
                list_name = atm_list[i].ele+'-'+atm_list[j].ele
                if list_name in dist_list:
                    dist_list[list_name].append(dist_matrix[i][j])
                else:
                    dist_list[list_name] = [dist_matrix[i][j]]
                f.write('{:.3f}\t'.format(dist_matrix[i][j]))
            f.write('\n')

if args.plot:
    pairs = args.pair
    colors = ['b', 'r', 'g', 'c', 'm', 'y']
    plt.figure(figsize=(10, 6))
    for idx, pair in enumerate(pairs):
        if pair in dist_list:
            distances = np.array(dist_list[pair])
            kde = gaussian_kde(distances, bw_method=args.bandwidth)
            distance_range = np.linspace(1.0, 7.0, 500)
            plt.plot(distance_range, kde(distance_range), label=f'{pair} bonds', color=colors[idx])
        else:
            raise ValueError('Pair is not in your POSCAR, please check.')
    plt.xlabel(r'Bond Length / $\AA$')
    plt.ylabel('Density / Arbitary Unit')
    plt.title('Pair Distribution')
    plt.grid(False)
    plt.legend()
    plt.xlim([1.0,7.0])
    plt.ylim([0.0,1.5])
    plt.savefig('bond_dist.png',transparent = True,format='png')
    plt.show()
        


