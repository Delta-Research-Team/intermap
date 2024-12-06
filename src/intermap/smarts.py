# Created by gonzalezroy at 12/6/24
import prolif.interactions as prolinter
from rdkit import Chem


def get_hydroph_smarts():
    """
    Get the hydrophobic SMARTS pattern

    Returns:
        smart_hydroph (str): SMARTS pattern of the hydrophobic interaction
    """
    lig_hydroph = prolinter.Hydrophobic().lig_pattern
    smart_hydroph = Chem.MolToSmarts(lig_hydroph)
    return smart_hydroph


def get_hb_smarts():
    """
    Get the hydrogen bond SMARTS pattern

    Returns:
        smart_acc (str): SMARTS pattern of the hydrogen bond acceptor interaction
        smart_don (str): SMARTS pattern of the hydrogen bond donor interaction
    """
    acceptor_hb = prolinter.HBAcceptor().lig_pattern
    donnor_hb = prolinter.HBDonor().prot_pattern
    smart_acc = Chem.MolToSmarts(acceptor_hb)
    smart_don = Chem.MolToSmarts(donnor_hb)
    return smart_acc, smart_don


def get_xb_smarts():
    """
    Get the halogen bond SMARTS pattern

    Returns:
        smart_acc (str): SMARTS pattern of the halogen bond acceptor interaction
        smart_don (str): SMARTS pattern of the halogen bond donor interaction
    """
    acceptor_xb = prolinter.XBAcceptor().lig_pattern
    donnor_xb = prolinter.XBAcceptor().prot_pattern
    smart_acc = Chem.MolToSmarts(acceptor_xb)
    smart_don = Chem.MolToSmarts(donnor_xb)
    return smart_acc, smart_don


def get_ionic_smarts():
    """
    Get the (cat/an)ionic SMARTS patterns

    Returns:
        smart_cation (str): SMARTS pattern of the cationic interaction
        smart_anion (str): SMARTS pattern of the anionic interaction
    """
    cations = prolinter.Cationic().lig_pattern
    anions = prolinter.Cationic().prot_pattern
    smart_cation = Chem.MolToSmarts(cations)
    smart_anion = Chem.MolToSmarts(anions)
    return smart_cation, smart_anion


def get_aromatic_rings():
    """
    Get the aromatic ring SMARTS pattern

    Returns:
        smart_rings (list): SMARTS pattern of the aromatic ring interaction
    """
    rings = prolinter.CationPi().pi_ring
    smart_rings = [Chem.MolToSmarts(ring) for ring in rings]
    return smart_rings


def get_metalic():
    """
    Get the metalic SMARTS pattern of the ligand

    Returns:
        smart_metal_acc (str): SMARTS pattern of the metalic acceptor interaction
        smart_metal_don (str): SMARTS pattern of the metalic donor interaction
    """
    metal_don = prolinter.MetalDonor().lig_pattern
    metal_acc = prolinter.MetalDonor().prot_pattern
    smart_metal_don = Chem.MolToSmarts(metal_don)
    smart_metal_acc = Chem.MolToSmarts(metal_acc)
    return smart_metal_acc, smart_metal_don


def get_first_residues(u):
    """
    Get the first residue of each type in the Universe

    Args:
        u (mda.Universe): Universe object

    Returns:
        first_occurences (dict): Dictionary with the first
         occurrence of each residue type
    """
    unique_resnames = set(u.residues.resnames)
    first_occurences = {}
    while unique_resnames:
        all_resids = u.residues
        for res in all_resids:
            if res.resname in unique_resnames:
                first_occurences[res.resname] = res
                unique_resnames.remove(res.resname)
                break
    return first_occurences


def query_mol(mol, query):
    """
    Get the atom names of the query in the molecule

    Args:
        mol (rdkit.Chem.rdchem.Mol): Molecule object
        query (rdkit.Chem.rdchem.Mol): Query object

    Returns:
        atom_names (str): Atom names of the query in the molecule
    """
    matches = mol.GetSubstructMatches(query)
    atom_names = ' '.join(res.atoms.names[[y for x in matches for y in x]])
    return atom_names


class ProlifSmarts:
    """
    Class to get the SMARTS patterns of the ligand interactions from prolif
    """

    def __init__(self):
        # Hydrophobics
        self.hydroph = get_hydroph_smarts()

        # Hydrogen bonds
        self.hb_acc, self.hb_don = get_hb_smarts()

        # Halogen bonds
        self.xb_acc, self.xb_don = get_xb_smarts()

        # Cations & Anions
        self.cations, self.anions = get_ionic_smarts()

        # Aromatic rings
        self.rings = get_aromatic_rings()

        # Metalic
        self.metal_acc, self.metal_don = get_metalic()

        # As a dictionary
        self.interactions = {
            "hydroph": self.hydroph,
            "hb_acc": self.hb_acc,
            "hb_don": self.hb_don,
            "xb_acc": self.xb_acc,
            "xb_don": self.xb_don,
            "cations": self.cations,
            "anions": self.anions,
            "rings6": self.rings[0],
            "rings5": self.rings[1],
            "metal_acc": self.metal_acc,
            "metal_don": self.metal_don
        }


# =============================================================================
# Debugging Smarts
# =============================================================================
self = ProlifSmarts()
smarts = self.interactions

# %% ==========================================================================
# Debugging sels
# =============================================================================
import MDAnalysis as mda

topo = '/media/gonzalezroy/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1.prmtop'
traj = '/media/gonzalezroy/Expansion/RoyData/oxo-8/raw/water/A2/8oxoGA2_1_sk100.nc'
u = mda.Universe(topo, traj)
u.select('protein')

unique_residues = get_first_residues(u)
t6 ='[a&r6]1:[a&r6]:[a&r6]:[a&r6]:[a&r6]:[a&r6]:1'
for case in unique_residues:
    if case == 'WAT':
        continue

    res = unique_residues[case]
    mol = res.atoms.convert_to("RDKIT", force=True)
    ri = mol.GetRingInfo()
    print(ri.AtomRings())

    query = Chem.MolFromSmarts(t6)
    matches = mol.GetSubstructMatches(query)
    print(case, matches)



    atom_names = query_mol(mol, query)
    selection = 'resname {} and name {}'.format(case, atom_names)



def isRingAromatic(mol, bondRing):
        for id in bondRing:
            if not mol.GetBondWithIdx(id).GetIsAromatic():
                return False
        return True
ri = mol.GetRingInfo()
print(isRingAromatic(mol, ri.BondRings()[2]))
