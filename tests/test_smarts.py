# Created by rglez at 2/2/25

# def get_singles(smart, universe, u_res, u_mol):
#     atom_names = u_res.atoms.names
#     query = Chem.MolFromSmarts(smart)
#     match = [list(x) for x in u_mol.GetSubstructMatches(query)]
#     match_names = [' '.join(atom_names[x]) for x in match]
#
#     singles = set()
#     for sel in match_names:
#         to_found = f'name {sel}'
#         idx = universe.select_atoms(to_found).indices
#         singles.update(idx)
#
#     return np.asarray(list(singles))
#
#
# topo1, trj1 = topo_trajs_fix['trj1']
#
# iman = IndexManager(topo1, trj1, 'all', 'all', 'all')
# universe = iman.universe
# renamed_universe = iman.renamed_universe
# u_mol = universe.atoms.convert_to("RDKIT", force=False)
# u_res = universe.residues
# unique_residues, unique_rdmols, unique_idx = get_uniques(universe)
#
# hp_imap = iman.get_singles('hydroph', iman.sel_idx)
# hp_mda = get_singles(smarts['hydroph'], universe, u_res, u_mol)
#
# np.setdiff1d(hp_mda, hp_imap)
