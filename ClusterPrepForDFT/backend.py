import numpy as np
import pymatgen as mg
from ClusterPrepForDFT import center_cluster_on_atom


class MolecularCluster(mg.core.Molecule):
    
    def __init__(self, species, coords, **kwargs):
        super().__init__(species, coords, **kwargs)
        
    @classmethod
    def from_file(cls, filename, cluster_radius, center_atom):
        structure = mg.core.Structure.from_file(filename)
        index = structure.species.index(mg.core.Element(center_atom))
        sites = structure.get_sites_in_sphere(structure[index].coords, cluster_radius)
        sites = [s[0] for s in sites] # drop calculated distances
        mol = cls.from_sites(sites)
        mol.translate_sites(vector=-structure[index].coords)
        return mol
        
    def display_jupyter(self, center_color='rgb(3,248,252)'):
        import nglview
        view = nglview.show_pymatgen(self)
        view.center()
        try:
            center_index = self.index(self.get_sites_in_sphere([0, 0, 0], 0.1)[0])
        except IndexError:
            center_index = None
        view.clear_representations()
        view.add_representation('ball+stick', selection='*')
        if center_color is not None:
            view.add_representation('ball+stick', selection=[center_index], color=center_color, radius=0.25)
        return view
