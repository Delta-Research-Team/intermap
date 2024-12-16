# Created by rglez at 12/9/24
"""
Cutoff distances for the different definitions of interactions
"""

# ==== Close contacts =========================================================
dist_cut_CloseContacts = 0.30

# ==== Ionics =================================================================
dist_cut_Ionic = 0.45

# ==== Hydrophobic ============================================================
dist_cut_Hydroph = 0.45

# ==== Metalic ================================================================
dist_cut_Metalic = 0.28

# ==== HBonds =================================================================
dist_cut_HA = 0.25
dist_cut_DA = 0.39
ang_cut_DHA = 90

# ==== XBonds =================================================================
dist_cut_XA = 0.25
dist_cut_XD = 0.39
ang_cut_DXA = 90

# ==== Aromatics ==============================================================
dist_cut_PiCation = 0.45
min_ang_PiCation = 0
max_ang_PiCation = 30

dist_cut_PiStacking = 0.6
min_dist_PiStacking = 0.38
min_ang_PiStacking = 0
max_ang_PiStacking = 90

dist_cut_EdgeToFace = 0.6
min_dist_EdgeToFace = 0.38
min_ang_EdgeToFace = 0
max_ang_EdgeToFace = 90

dist_cut_FaceToFace = 0.45
min_dist_FaceToFace = 0.38
min_ang_FaceToFace = 0
max_ang_FaceToFace = 90
