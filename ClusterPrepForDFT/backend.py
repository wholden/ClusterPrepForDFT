import numpy as np
import pymatgen as mg
from ClusterPrepForDFT import center_cluster_on_atom


class MolecularCluster(mg.core.Molecule):
    
    def __init__(self, species, coords, **kwargs):
        super().__init__(species, coords, **kwargs)

    @staticmethod
    def process_structure_to_molecular_cluster(structure, cluster_radius, center_atom):
        index = structure.species.index(mg.core.Element(center_atom))
        sites = structure.get_sites_in_sphere(structure[index].coords, cluster_radius)
        sites = [s[0] for s in sites] # drop calculated distances
        new = MolecularCluster.from_sites(sites)
        new.translate_sites(vector=-structure[index].coords)
        return new
        
    @classmethod
    def from_file(cls, filename, cluster_radius, center_atom):
        structure = mg.core.Structure.from_file(filename)
        new = cls.process_structure_to_molecular_cluster(structure, cluster_radius, center_atom)
        new._structure = structure
        return new

    # @classmethod
    # def from_dict(cls, d, cluster_radius, center_atom):
    #     structure = mg.core.Structure.from_dict(d)
    #     new = cls.process_structure_to_molecular_cluster(structure, cluster_radius, center_atom)
    #     return new

    # @classmethod
    # def from_sites(cls, sites, cluster_radius=0, center_atom='', skip_process=False):
    #     if not skip_process:
    #         structure = mg.core.Structure.from_sites(sites)
    #         new = cls.process_structure_to_molecular_cluster(structure, cluster_radius, center_atom)
    #     else:
    #         new = super().from_sites(sites)
    #     return new

    # @classmethod
    # def from_str(cls, string, cluster_radius, center_atom):
    #     structure = mg.core.Structure.from_str(string)
    #     new = cls.process_structure_to_molecular_cluster(structure, cluster_radius, center_atom)
    #     return new

    @staticmethod
    def classify_unique_sites(structure):
        uniques = []
        for site in structure:
            parsed_props = parse_site_props(structure, structure.index(site))
            site._parsed_props = parsed_props
            if parsed_props not in uniques:
                uniques.append(parsed_props)
        for site in structure:
            site.classification = uniques.index(site._parsed_props)
        return len(uniques), uniques

    @staticmethod
    def parse_site_props(structure, index):
        cnn = CrystalNN()
        props = {}
        try:
            props['element'] = str(structure[index].specie.element)
            props['oxid-state'] = structure[index].specie.oxi_state
        except AttributeError:
            props['element'] = str(structure[index].element)
            props['oxid-state'] = None
        neighbors = cnn.get_nn(structure, index)
        neighbors = [n.species_string for n in neighbors]
        props['neighbors'] = {s:neighbors.count(s) for s in set(neighbors)}
        return props
        
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
