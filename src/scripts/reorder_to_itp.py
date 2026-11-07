"""Write an SDF of SMILES in the atom order of a LigParGen ITP (graph isomorphism)."""
import sys, re
import networkx as nx
from networkx.algorithms import isomorphism
from rdkit import Chem
from rdkit.Chem import AllChem
smiles, itp, out = sys.argv[1:4]
text = open(itp).read()
atoms = re.search(r"\[ atoms \](.*?)\[ bonds \]", text, re.S).group(1)
elem = {}
for line in atoms.splitlines():
    f = line.split()
    if len(f) >= 8 and f[0].isdigit():
        elem[int(f[0])] = 'H' if float(f[7]) < 2 else {12: 'C', 14: 'N', 15: 'O', 16: 'O', 32: 'S', 35: 'Cl', 19: 'F'}[round(float(f[7]))]
bonds = re.search(r"\[ bonds \](.*?)\[\s*angles", text, re.S).group(1)
G = nx.Graph()
for i, e in elem.items(): G.add_node(i, e=e)
for line in bonds.splitlines():
    f = line.split()
    if len(f) >= 2 and f[0].isdigit(): G.add_edge(int(f[0]), int(f[1]))
mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
H = nx.Graph()
for a in mol.GetAtoms(): H.add_node(a.GetIdx(), e=a.GetSymbol())
for b in mol.GetBonds(): H.add_edge(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
gm = isomorphism.GraphMatcher(G, H, node_match=lambda a, b: a['e'] == b['e'])
mapping = next(gm.isomorphisms_iter())  # itp index (1-based) -> rdkit index
order = [mapping[i] for i in sorted(elem)]
mol = Chem.RenumberAtoms(mol, order)
params = AllChem.ETKDGv3(); params.randomSeed = 1
assert AllChem.EmbedMolecule(mol, params) == 0
AllChem.MMFFOptimizeMolecule(mol)
Chem.MolToMolFile(mol, out)
print("wrote", out)
