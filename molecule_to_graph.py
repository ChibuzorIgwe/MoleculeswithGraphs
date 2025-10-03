from rdkit import Chem
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import AllChem
import py3Dmol
import cirpy
import torch
from torch_geometric.utils import from_smiles
import tempfile
import webbrowser
import os


def mol2graph(molecule):
   
    molecule_smile = cirpy.resolve(molecule, "smiles") 
    data = from_smiles(molecule_smile, with_hydrogen=True) #this extracts the graph representation of the molecule
    print(f"Then name of the molecule is: {molecule}")
    print(f"The summary of molecule's 'data' is: {data}") 
    print(f"The number of nodes in graph is: {data.x.shape[0]}") 
    print(f"The number of edges in graph is: {data.edge_index.t().shape[0]}") 
    print(f"The number of node features in graph is: {data.x.shape[1]}") 
    print(f"The number of edge features in graph is: {data.edge_attr.shape[1]}") 
    print(f"The shape of the node feature matrix is: {data.x.shape}") 
    print(f"The shape of the edge index list/ or two-column matrix is: {data.edge_index.t().shape}") 
    print(f"The shape of the edge attribute matrix is: {data.edge_attr.shape}") 
    print(f"The node features of an atom in the graph are: {data.x[0]}")  
    print(f"The node features of the graph are: {data.x}") 

    mol = Chem.MolFromSmiles(molecule_smile)
    mol

    atom = mol.GetAtomWithIdx(0) # Get the first atom (a Carbon)

    print(f"Atomic Number: {atom.GetAtomicNum()}")
    print(f"Atom Index: {atom.GetIdx()}")
    print(f"Atom Mass: {atom.GetMass()}")
    print(f"Atom Degree: {atom.GetDegree()}")
    print(f"Atom Hybridization: {atom.GetHybridization()}")
    print(f"Atom Is Aromatic: {atom.GetIsAromatic()}")
    print(f"Atom Is In Ring: {atom.IsInRing()}")
    print(f"Atom Type: {atom.GetSymbol()}")
    print(f"Formal Charge: {atom.GetFormalCharge()}")
    print(f"Chirality: {atom.GetChiralTag()}")
    print(f"Number of Hydrogens: {atom.GetTotalNumHs()}")

    mol3d = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol3d, AllChem.ETKDGv3())
    AllChem.UFFOptimizeMolecule(mol3d)
    mb = Chem.MolToMolBlock(mol3d)

    # Create HTML content for 3D visualization since py3Dmol doesn't work unless its in an interactive python notebook
    html_content = f"""<!DOCTYPE html>

<html>
<head>
    <title>3D Molecule Viewer</title>
    <script src="https://cdn.jsdelivr.net/npm/3dmol@2.5.2/build/3Dmol-min.js"></script>
</head>
<body>
<h3>{molecule}</h3>
    <div id="viewer" style="width: 400px; height: 300px;"></div>
    <script>
        let viewer = $3Dmol.createViewer(document.getElementById('viewer'), {{
            backgroundColor: 'white'
        }});
        viewer.addModel(`{mb}`, "sdf");
        viewer.setStyle({{ elem: 'C' }}, {{
            stick: {{ color: 'red', radius: 0.2 }},
            sphere: {{ color: 'red', scale: 0.35 }}
        }});
        viewer.setStyle({{ elem: 'H' }}, {{
            stick: {{ color: 'white', radius: 0.2 }},
            sphere: {{ color: 'white', scale: 0.35 }}
        }});
        viewer.addLabel('{molecule}', {{ position: {{ x: 3, y:3, z: 3 }}, backgroundColor: 'rgba(255,255,255,0.8)', fontColor: 'black', fontSize: 24, showBackground: true, inFront: true}});
        viewer.zoomTo();
        viewer.render();
    </script>
</body>
</html>"""
        
    # Write to temporary file
    with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as f:
        f.write(html_content.encode("utf-8"))
        temp_path = f.name
        
    # Open in default browser
    webbrowser.open(f"file://{os.path.abspath(temp_path)}")

if __name__ == "__main__":

    # List of alkanes from C1 to C10
    molecules = ["methane",  "ethane",  "propane",  "butane",  "pentane",  "hexane",  "heptane",  "octane",  "nonane",  "decane"]

    # Visualize each alkane molecule
    for molecule in molecules:
        mol2graph(molecule)
