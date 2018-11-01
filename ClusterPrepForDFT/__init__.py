import numpy as np
import os
import pymatgen as mg


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
            bondcharges[k+k1] = float(k[-1]+k[-2]) / v1 # no 10+
            
    totalcharge = 0
    for k, v in missingbonds.items():
        totalcharge += v * bondcharges[k]
    return totalcharge, bondcharges


def center_cluster_on_atom(Molecule, atom):
    r = 0.1
    centeratom = atom
    while True:
        sites = Molecule.get_sites_in_sphere(Molecule.center_of_mass, r)
        found = [s[0].specie.element.name == centeratom for s in sites]
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
