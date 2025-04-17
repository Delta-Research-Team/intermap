"""
Helper functions for the InterMap Visualizations app.

"""

import base64
import pandas as pd
from pathlib import Path


def find_topology_file(temp_path, original_path):

    TOPOLOGY_EXTENSIONS = {
        '.pdb', '.gro', '.prmtop', '.top', '.psf', '.mol2',
        '.xyz', '.pqr', '.gsd', '.data', '.lammpstrj'
    }

    try:
        original_dir = Path(original_path).parent
        base_name = Path(original_path).stem

        print(
            f"Searching for topology file in original directory: {original_dir}")
        print(f"Base name of the original file: {base_name}")

        for ext in TOPOLOGY_EXTENSIONS:
            potential_topo = original_dir / f"{base_name}{ext}"
            if potential_topo.exists():
                print(
                    f"Topology file found: {potential_topo}")
                return True, str(potential_topo)

        for file in original_dir.iterdir():
            if file.suffix.lower() in TOPOLOGY_EXTENSIONS:
                print(
                    f"Topology file found: {file}")
                return True, str(file)


        print("No topology file found.")
        return False, ""

    except Exception as e:
        print(f"Error searching for topology file: {str(e)}")
        return False, ""

def get_image_base64(image_path):
    """
    Convert an image file to base64 string.

    Args:
        image_path (str): Path to the image file

    Returns:
        str: Base64 encoded image string
    """
    try:
        with open(image_path, "rb") as image_file:
            encoded_string = base64.b64encode(image_file.read()).decode()
            return f"data:image/jpeg;base64,{encoded_string}"
    except Exception as e:
        print(f"Error loading image: {e}")
        return ""

def calculate_prevalence(row):
    """
    Calculate the prevalence of an interaction.

    Args:
        row (pd.Series): Row containing frame data

    Returns:
        float: Percentage of frames where the interaction occurs
    """
    if 'prevalence' in row.index:
        return float(row['prevalence'])

def generate_interaction_choices(df):
    """
    Generate choices for interaction type checkboxes.

    Args:
        df (pd.DataFrame): DataFrame containing interaction data

    Returns:
        dict: Dictionary of interaction types with their counts
    """
    if df is None:
        return {}

    interaction_counts = df['interaction_name'].value_counts()
    choices = {
        interaction: f"{interaction} ({count})"
        for interaction, count in interaction_counts.items()
    }
    return dict(sorted(choices.items()))

def validate_mda_selection(selection):

    VALID_KEYWORDS = {
        'resname', 'name', 'resid', 'backbone', 'protein', 'water',
        'nucleic', 'segment', 'segid', 'chain', 'residue', 'around',
        'within', 'same', 'byres', 'all'
    }

    if not selection.strip():
        return False, "The selection is empty."

    invalid_chars = set('!@#$%^&*+=[]{}|\\;:"\'<>?')
    if any(char in invalid_chars for char in selection):
        return False, "The selection contains invalid characters."

    words = set(selection.lower().split())
    if not any(keyword in words for keyword in VALID_KEYWORDS):
        return False, "The selection does not contain valid MDAnalysis keywords."

    return True, ""

