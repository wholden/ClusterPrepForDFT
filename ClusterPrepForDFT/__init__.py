import numpy as np
import os
import pymatgen as mg
import itertools
from scipy.optimize import minimize


def make_molecular_cluster_from_cif(file, supercell=None, guessoxid=False):
    structure = mg.Structure.from_file(file)
    if guessoxid:
        structure.add_oxidation_state_by_guess()
    if supercell is not None:
        structure.make_supercell(supercell)
        
    molecule = mg.Molecule(structure.species, structure.cart_coords)
    return molecule


def assign_nearest_neighbors(Molecule, rscale=1.1, maxr=5):
    rscale = 1.1
    for m in Molecule:
        neighbors = Molecule.get_neighbors(m, 5)
        dists = [n[1] for n in neighbors]
        withinr = np.min(dists) * rscale

        nearests = []
        for n, d in neighbors:
            if d < withinr:
                nearests.append(n)
        m.nearests = nearests


def calculate_missing_bonds(Molecule, bondmap):
    totalmissing = {}
    for m in Molecule:
        missing = bondmap[m.species_string].copy()
        for n in m.nearests:
            if n.species_string in missing:
                missing[n.species_string] -= 1
        for k, v in missing.items():
            try:
                totalmissing[m.species_string+k] += v
            except KeyError:
                totalmissing[m.species_string+k] = v
    return totalmissing


def calculate_charge_for_cluster(missingbonds, bondmap):
    bondcharges = {}
    for k, v in bondmap.items():
        for k1, v1 in v.items():
            try:
                bondcharges[k+k1] = float(k[-1]+k[-2]) / v1 # no 10+
            except ValueError: # handle single minus case
                print('Charge parse error, trying case +/- 1')
                bondcharges[k+k1] = float(k[-1]+'1') / v1 # no 10+
            
    totalcharge = 0
    for k, v in missingbonds.items():
        totalcharge += v * bondcharges[k]
    return totalcharge, bondcharges


def center_cluster_on_atom(Molecule, atom):
    r = 0.1
    centeratom = atom
    while True:
        sites = Molecule.get_sites_in_sphere(Molecule.center_of_mass, r)
        try:
            found = [s[0].specie.element.name == centeratom for s in sites]
        except AttributeError:
            found = [s[0].specie.name == centeratom for s in sites]
        if any(found):
            break
        r += 0.1

    centerindex = np.where(found)[0][0]
    center = sites[centerindex][0]
    Molecule.translate_sites(vector=-center.coords)


def dumb_rename_center_atom(file, Molecule):
    index = np.where([all(el == np.array([0,0,0])) for el in Molecule.cart_coords])[0][0] + 2
    with open(file, 'r+') as f:
        lines = f.readlines()
        lines[index] = '{}1 {} {} {}\n'.format(*lines[index].split())
        f.seek(0)
        f.writelines(lines)


def cif_to_charged_cluster(file, bondmap, centeratom=False, supercell=None, oxid='guess', clobber=False):
    assert file[-4:] == '.cif', 'Only works with cif files'
    if oxid == 'guess':
        mol = make_molecular_cluster_from_cif(file, supercell=supercell, guessoxid=True)
    else:
        mol = make_molecular_cluster_from_cif(file, supercell=supercell, guessoxid=False)
        mol.add_oxidation_state_by_element(oxid)
    
    assign_nearest_neighbors(mol)
    missing = calculate_missing_bonds(mol, bondmap)
    charge, bondcharges = calculate_charge_for_cluster(missing, bondmap)
    
    if centeratom:
        center_cluster_on_atom(mol, centeratom)
    
    outfilename = file[:-4]+'.xyz'
    
    if not clobber:
        assert not os.path.exists(outfilename), 'Output file with name {} already exists!'.format(outfilename)

    mol.add_oxidation_state_by_element(dict.fromkeys(mol.symbol_set))
    mol.to('xyz', outfilename)

    with open(outfilename, 'r+') as file:
        lines = file.readlines()
        lines[1] = '{}; charge for DFT cluster calculation: {}; based on bondchargemap: {};\n'.format(lines[1][:-1], charge, bondcharges)
        file.seek(0)
        file.writelines(lines)

    dumb_rename_center_atom(outfilename, mol)


def cif_to_charged_cluster_within_radius(file, bondmap, radius, centeratom=False, supercell=None, oxid='guess', clobber=False):
    '''Heads up, supercell needs to be large enough to cover the radius requested'''
    assert file[-4:] == '.cif', 'Only works with cif files'
    if oxid == 'guess':
        mol = make_molecular_cluster_from_cif(file, supercell=supercell, guessoxid=True)
    else:
        mol = make_molecular_cluster_from_cif(file, supercell=supercell, guessoxid=False)
        mol.add_oxidation_state_by_element(oxid)
    
    center_cluster_on_atom(mol, centeratom)
    sites = mol.get_sites_in_sphere([0, 0, 0], radius)
    mol = mg.Molecule([s[0].species_string for s in sites], [s[0].coords for s in sites])
    assign_nearest_neighbors(mol)
    
    missing = calculate_missing_bonds(mol, bondmap)
    charge, bondcharges = calculate_charge_for_cluster(missing, bondmap)
    
        
    
    outfilename = file[:-4]+'.xyz'
    
    if not clobber:
        assert not os.path.exists(outfilename), 'Output file with name {} already exists!'.format(outfilename)

    mol.add_oxidation_state_by_element(dict.fromkeys(mol.symbol_set))
    mol.to('xyz', outfilename)

    with open(outfilename, 'r+') as file:
        lines = file.readlines()
        lines[1] = '{}; charge for DFT cluster calculation: {}; based on bondchargemap: {};\n'.format(lines[1][:-1], charge, bondcharges)
        file.seek(0)
        file.writelines(lines)

    dumb_rename_center_atom(outfilename, mol)


def make_bondmap_from_cif():
    #TODO E
    pass


#Example usage
# m = cif_to_charged_cluster(
#     'V2O3_PDFCard-01-071-0280.cif', 
#     {'V3+':{'O2-': 6},'O2-':{'V3+': 4}}, 
#     centeratom='V',
#     supercell=1, 
#     clobber=True)

def del_atom_at_coord(Molecule, coord):
	#Example usage
	# del_atom_at_coord(mol, [-5.68529,-3.53916,11.34309])
    sites = Molecule.get_sites_in_sphere(coord, 0.001)
    assert len(sites) == 1, 'Multiple sites found at coords'
    site = sites[0][0]
    siteindex = Molecule.index(site)
    print('Removing site {} at index {}'.format(site, siteindex))
    Molecule.remove_sites([siteindex])


def get_nearest_sites_by_num_neighbors(mol, site, numnearest):
    num = numnearest[str(site.specie)]
    
    neighbors = []
    rinit= 0.4 # angstroms
    while len(neighbors) < num:
        rinit += 0.1
        neighbors = mol.get_neighbors(site, rinit)
        if rinit > 10:
            raise RuntimeError('Number of requested nearest neighbors not found')
    return neighbors


def make_hydrogen_capped_cluster_from_cif(outfilename, ciffile, numnearest, cappingdistances, bqchargemap, centeratom, sphereradius=6, supercell=6):
    '''Make hydrogen-capped cluster.
    
    Supercell needs to be big enough that all atoms in sphere and neighbors thereof
    have full neighbors.'''
    print('Importing structure')
    struct = mg.Structure.from_file(ciffile)
    struct.make_supercell(6)
    mol = mg.Molecule(struct.species, struct.cart_coords)
    center_cluster_on_atom(mol, centeratom) # fixed 3.23.19 to add centeratom arg

    print('Finding site inside/outside sphere.')
    innersites = mol.get_sites_in_sphere([0, 0, 0], sphereradius)
    inner = mg.Molecule.from_sites([s[0] for s in innersites])
    SHELLDIST = 4
    outersites = mol.get_neighbors_in_shell([0, 0, 0], sphereradius + SHELLDIST / 2, SHELLDIST / 2)

    outersiteindices = []
    for site, _ in outersites:
        outersiteindices.append(mol.index(site))
    neighboringsites = []
    for osi in outersiteindices:
        print('Checking neighbors of site {}'.format(osi), end='\r')
        neighbors = get_nearest_sites_by_num_neighbors(mol, mol[osi], numnearest)
        for s, _ in neighbors:
            if s in inner:            
                try:
                    mol[osi].innerneighbors.append(s)
                except AttributeError:
                    mol[osi].innerneighbors = [s]

                if osi not in neighboringsites:
                    neighboringsites.append(osi)

    hydrcoords = []
    hydrneighbors = []
    hydrreplace = []
    for ns in neighboringsites:
        ocoords = mol[ns].coords
        for n in mol[ns].innerneighbors:
            hydrcoords.append(np.mean([ocoords, n.coords], axis=0))
            hydrneighbors.append(n)
            hydrreplace.append(mol[ns])

    hydrsites = [mg.Site('H', c) for c in hydrcoords]
    hydr = mg.Molecule.from_sites(hydrsites)

    # add hydrogen neighbors and adjust bond lengths accordingly
    for hs, hn, hr in zip(hydrsites, hydrneighbors, hydrreplace):
        hydr[hydr.index(hs)].innerneighbor = hn
        hydr[hydr.index(hs)].replacing = hr
    for s in hydr:
        adjust_bond_length(s, s.innerneighbor, cappingdistances[str(s.innerneighbor.specie)])

    capped = mg.Molecule.from_sites(hydr.sites + inner.sites)
    print('Writing capped cluster file: ', '{}.xyz'.format(outfilename))
    capped.to('xyz', '{}.xyz'.format(outfilename))
    print('Writing bq charge file: ', '{}.bq'.format(outfilename))
    with open('{}.bq'.format(outfilename), 'w') as file:
        for s in hydr:
            file.write('Bq {:>15.8f}{:>15.8f}{:>15.8f}{:>10.4f}\n'.format(
                *s.coords, bqchargemap[str(s.replacing.specie)]))


def make_constrain_ordered_xyz(file, neighbor0='P', neighbor1='H'):
    mol = mg.Molecule.from_file(file)
    assign_nearest_neighbors(mol)

    edge = mg.Molecule(['Zr'], [[0,0,0]]) #dummy Zr to let me instantiate the molecule
    edge.pop(0) #discard dummy

    poplist = []
    for i, site in enumerate(mol):
        if site.specie.name == neighbor0:
            for n in site.nearests:
                if n.specie.name == neighbor1:
                    poplist.append(i)
                    edge.append(site.specie, site.coords)
                    continue
    mol.remove_sites(poplist)

    poplist = []
    for i, site in enumerate(mol):
        if site.specie.name == neighbor1:
            for n in site.nearests:
                if n.specie.name == neighbor0:
                    poplist.append(i)
                    edge.append(site.specie, site.coords)
                    continue
    mol.remove_sites(poplist)

    filename = file.split('.xyz')[0].split('/')[-1] #oops only works on linux...
#     unconstrainf = '{}_Unconstrain.xyz'.format(filename)
#     mol.to('xyz', unconstrainf)
#     constrainf = '{}_Constrain.xyz'.format(filename)
#     edge.to('xyz', constrainf)

#     unca, uncc = read_xyz(unconstrainf)
#     ca, cc = read_xyz(constrainf)

    unca, uncc = mol.species, mol.cart_coords
    ca, cc = edge.species, edge.cart_coords

    print('Saving new file as: {}'.format('{}_Reordered'.format(filename)))
    write_xyz_from_atoms_coords('{}_Reordered'.format(filename), 
                                atoms=unca + ca, 
                                coords=np.concatenate([uncc, cc]), 
                                comment='Unconstrained {}; Constrained {}; Fix atoms {}:{}'.format(
                                    len(unca),
                                    len(ca),
                                    len(unca) + 1,
                                    len(unca) + len(ca)
                                ))


def get_neighbor_site_angles(mol, site, r=2.5):
    neighbors = [n[0] for n in mol.get_neighbors(site, r)]
    angles = []
    for comb in itertools.combinations(neighbors, 2):
        j = mol.index(site)
        i, k = [mol.index(s) for s in comb]
        angles.append(mol.get_angle(i, j, k))
    return sorted(np.array(angles))


def adjust_bond_length(site, neighbor, targetdist):
    directionvec = site.coords - neighbor.coords
    directionvec /= np.linalg.norm(directionvec)
    oldcoords = site.coords
    def func(c):
        site._coords = oldcoords + c * directionvec
        return (site.distance(neighbor) - targetdist) ** 2
    minimize(func, 0.1)


def min_sum_squared(vec0, vec1):
    sumsqs = []
    for p in itertools.permutations(vec1):
        sumsqs.append(sum((vec0 - np.array(p))**2))
    return np.min(sumsqs)


def avg_sum_squared(vec0, vec1):
    sumsqs = []
    for p in itertools.permutations(vec1):
        sumsqs.append(sum((vec0 - np.array(p))**2))
    return np.mean(sumsqs)


def total_sum_squared(vec0, vec1):
    sumsqs = []
    for p in itertools.permutations(vec1):
        sumsqs.append(sum((vec0 - np.array(p))**2))
    return np.sum(sumsqs)


def adjust_site_to_target_angles(mol, adjustsite, neighborsite, targetangles, r=2.5, minfunc=avg_sum_squared):
    def func(coords):
        adjustsite._coords = coords
        return minfunc(targetangles, get_neighbor_site_angles(mol, neighborsite, r))
    minimize(func, adjustsite._coords)
    return adjustsite._coords
