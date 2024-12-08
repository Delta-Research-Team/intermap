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

            # single atoms matches
            "hydroph": self.hydroph,
            "cations": self.cations,
            "anions": self.anions,
            "metal_acc": self.metal_acc,
            "metal_don": self.metal_don,
            "hb_acc": self.hb_acc,
            "xb_acc": self.xb_acc,

            # multiple atoms matches
            "hb_don": self.hb_don,
            "xb_don": self.xb_don,
            "rings5": '[a&r]1:[a&r]:[a&r]:[a&r]:[a&r]:1',
            "rings6": '[a&r]1:[a&r]:[a&r]:[a&r]:[a&r]:[a&r]:1'
        }

# =============================================================================
# Debugging Smarts
# =============================================================================
# self = ProlifSmarts()
# smarts = self.interactions
