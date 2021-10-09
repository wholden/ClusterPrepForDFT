import numpy as np
import pymatgen as mg
from ClusterPrepForDFT import center_cluster_on_atom


class MolecularCluster(mg.core.Molecule):
    
    def __init__(self, species, coords, **kwargs):
        super().__init__(species, coords, **kwargs)
        
    @classmethod
    def from_file(cls, filename, cluster_radius, center_atom):
        structure = mg.core.Structure.from_file(filename)
        cart_sizes = [
            structure.cart_coords[:, 0].max() - structure.cart_coords[:, 0].min(),
            structure.cart_coords[:, 1].max() - structure.cart_coords[:, 1].min(),
            structure.cart_coords[:, 2].max() - structure.cart_coords[:, 2].min(),
        ]
        structure.make_supercell([np.ceil(cluster_radius * 2 / sz) for sz in cart_sizes])
        mol = super().from_sites(structure.sites)
        center_cluster_on_atom(mol, center_atom)
        sites = mol.get_sites_in_sphere([0, 0, 0], cluster_radius)
        return super().from_sites(sites)
        
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
