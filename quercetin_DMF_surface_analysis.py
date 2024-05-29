#
# This script can be used for any purpose without limitation subject to the
# conditions at http://www.ccdc.cam.ac.uk/Community/Pages/Licences/v2.aspx
#
# This permission notice and the following statement of attribution must be
# included in all copies or substantial portions of this script.
#
# 2021-11-15: created by Alex Moldovan, The Cambridge Crystallographic Data Centre
# 2024-03-05: modified by Andy Maloney,  The Cambridge Crystallographic Data Centre
#

from ccdc import io, particle
import itertools
import numpy as np


def calculate_length(origin, target):
    if not isinstance(origin, np.ndarray) or not (isinstance(target, np.ndarray)):
        raise TypeError("Please supply numpy arrays for lengths.")
    return np.linalg.norm(target - origin)


def compute_triangle_area(a, b, c):
    s = (a + b + c) / 2
    return np.sqrt(s * (s - a) * (s - b) * (s - c))


class SurfaceChemistryAnalyser:
    """
    Calculates surface areas of different chemical moieties based on identification of atom types.
    """
    def __init__(self, surface_contacts, surface_chemistry, triangles):
        self.surface_contacts = surface_contacts
        self._chemistry = surface_chemistry
        self.triangles = triangles
        self.node_chemistry = None
        self.triangle_lookup = None
        self.chemistry_triangles = None
        self.node_search_dictionary = None

    def find_triangle_chemistries(self):
        """Calculates the chemistry near every triangle"""
        self.node_chemistry = self._get_node_chemistry(self.surface_contacts, self._chemistry)
        self.node_search_dictionary = self._make_node_search_dictionary(self.node_chemistry)
        self.triangle_lookup = self._triangle_lookup(self.node_search_dictionary, self.triangles)
        self._calculate_chemistry_lookup()

    def _calculate_chemistry_lookup(self):
        """Calculates a lookup table for the different chemistries with the triangle index"""
        if self.triangle_lookup is None:
            raise LookupError(f"Triangle look up not calculated. Wrong order of operations.")
        self._unique_triangle_lookup = self.unique_triangle_chemistries(self.triangle_lookup)
        self.chemistry_triangles = self._invert_list_to_dictionary(self._unique_triangle_lookup)

    def report_duplicated_chemistry(self, exclusion_list=None):
        """Returns atoms that occur in multiple fragments"""
        if not exclusion_list:
            exclusion_list = []
        if not isinstance(exclusion_list, list):
            raise ValueError("Supplied exclusion list is not in a list format.")
        return {atom: chemistry for atom, chemistry in
                self._invert_chemistry_dictionary(self._chemistry, exclusion_list).items() if len(chemistry) >= 2}

    @property
    def triangle_areas(self):
        """Calculates all the areas of a triangle"""
        return np.array(self.calculate_areas_of_triangles(self.triangles))

    @staticmethod
    def _get_node_chemistry(surface_contacts, chemistry):
        """Cross-references defined chemistry for each atom"""
        node_chemistry = {}
        for chem, atoms in chemistry.items():
            if len(atoms) == 0:
                continue
            for atom in atoms:
                if atom[-1] not in surface_contacts:
                    nodes = []
                else:
                    nodes = list(surface_contacts[atom[-1]])
                if chem not in node_chemistry.keys():
                    node_chemistry[chem] = nodes
                else:
                    node_chemistry[chem].extend(nodes)
        return node_chemistry

    @staticmethod
    def _make_node_search_dictionary(node_chemistry):
        """Makes a node search dictionary"""
        search_dictionary = {}
        for chem, values in node_chemistry.items():
            for v in values:
                if v not in search_dictionary.keys():
                    search_dictionary[v] = [chem]
                elif chem not in search_dictionary[v]:
                    search_dictionary[v].append(chem)
        return search_dictionary

    @staticmethod
    def _triangle_lookup(search_dictionary, triangles):
        """ Creates a list from the triangulated topology with each index containing the node chemistries"""
        new_triangles = []
        count = 0
        for triangle in triangles:
            new_triangle = []
            for n in triangle:
                try:
                    new_triangle.append(search_dictionary[n])
                except KeyError:
                    count += 1
            new_triangles.append(new_triangle)
        return new_triangles

    @staticmethod
    def unique_triangle_chemistries(triangles_lookup):
        """Stores unique data for triangle list with chemistries"""
        return [set(itertools.chain(*triangles)) for triangles in triangles_lookup]

    def area_of_chemistry(self, triangles_area, chemistry):
        """Returns area covered by specific chemistry in A^2"""
        try:
            chemistry_triangles = self.chemistry_triangles[chemistry]
        except KeyError:
            raise KeyError(
                f"Chemistry given is not listed in the look-up table. Select from the following types: "
                f"\n {self.chemistry_triangles.keys()}")
        return triangles_area[chemistry_triangles].sum()

    def percentage_coverage(self, triangles_area, chemistry):
        """Returns % coverage of a given chemistry from triangle area."""
        return (self.area_of_chemistry(triangles_area, chemistry) / triangles_area.sum()) * 100

    @staticmethod
    def _invert_chemistry_dictionary(d, exclusion_list=None):
        """Inverts the chemistry dictionary"""
        if not exclusion_list:
            exclusion_list = []
        inv_dict = {}
        for chem, atoms in d.items():
            if chem in exclusion_list:
                continue
            for atom in atoms:
                inv_dict.setdefault(atom[-1], set()).add(chem)
        return inv_dict

    @staticmethod
    def _invert_list_to_dictionary(n):
        """Inverts a list with sets to a dictionary"""
        new_dic = {}
        for i, values in enumerate(n):
            for j in values:
                new_dic.setdefault(j, []).append(i)
        return new_dic

    @staticmethod
    def calculate_areas_of_triangles(triangles):
        """Calculates area of individual triangles from node positions using Heron's formula"""
        triangle_areas = []
        for triangle in triangles:
            pos_0, pos_1, pos_2 = np.array(triangle[0]), np.array(triangle[1]), np.array(triangle[2])
            a_dist = calculate_length(pos_0, pos_1)
            b_dist = calculate_length(pos_0, pos_2)
            c_dist = calculate_length(pos_1, pos_2)
            triangle_areas.append(compute_triangle_area(a_dist, b_dist, c_dist))
        return triangle_areas

    @staticmethod
    def find_atom_types(surface, chemistry):
        """Returns the Sybyl Type for a given set of atoms."""
        groups = {sybyl_type: [] for sybyl_type in chemistry}
        for atom in surface.surface_atoms:
            groups[atom.sybyl_type].extend([[atom.coordinates, atom]])
        return groups


def analyse(input_crystal, hkl, offset):
    molecule = input_crystal.molecule
    molecule.assign_bond_types(which="unknown")
    input_crystal.molecule = molecule

    surface = particle.Surface(input_crystal, hkl, offset=offset)

    surface_chemistry = list(set([atom.sybyl_type for atom in surface.surface_atoms]))

    atom_types = SurfaceChemistryAnalyser.find_atom_types(surface, surface_chemistry)

    atom_surface_nodes_contact = surface.atom_surface_node_contacts

    atom_types_analysis = SurfaceChemistryAnalyser(atom_surface_nodes_contact,
                                                   atom_types,
                                                   list(surface.topology.triangles))

    atom_types_analysis.find_triangle_chemistries()
    # calculate the area occupied by each triangle element of the surface
    areas = atom_types_analysis.triangle_areas

    for atom_sybyl in atom_types.keys():
        print(f"{atom_sybyl} : {round(atom_types_analysis.percentage_coverage(areas, atom_sybyl), 3)}")
    analysis = {atom_sybyl: round(atom_types_analysis.percentage_coverage(areas, atom_sybyl), 3) for atom_sybyl in
                atom_types.keys()}

    return analysis


if __name__ == '__main__':
    # simple demonstration of how to use this script
    csd = io.CrystalReader('CSD')
    test_crystal = csd.crystal("HXACAN")
    test_crystal.assign_bonds()
    output = analyse(input_crystal=test_crystal, hkl=(0, 0, 2), offset=0.00)
    print(output)
