from rdkit import Chem
from rdkit.Chem import AllChem
import py3Dmol
import cirpy
import tempfile
import webbrowser
import os

def ballnstick_viz(molecule):
    """
    Visualize the molecule from a common name in 3D.
    Parameters:
        molecule (str): Common name of the molecule
    Returns:
        None: Displays the molecule in 3D
    """
    molecule_smile = cirpy.resolve(molecule, "smiles")
    mol = Chem.MolFromSmiles(molecule_smile)
    
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
    <div id="view" style="width: 450px; height: 350px;"></div>
    <script>
        let view = $3Dmol.createViewer(document.getElementById('view'), {{
            backgroundColor: 'white'
        }});
        view.addModel(`{mb}`, "sdf");
        view.setStyle({{ elem: 'C' }}, {{
            stick: {{ color: 'red', radius: 0.2 }},
            sphere: {{ color: 'red', scale: 0.35 }}
        }});
        view.setStyle({{ elem: 'H' }}, {{
            stick: {{ color: 'white', radius: 0.2 }},
            sphere: {{ color: 'white', scale: 0.35 }}
        }});
        view.addLabel('{molecule}', {{ position: {{ x: 3, y:3, z: 3 }}, backgroundColor: 'rgba(255,255,255,0.8)', fontColor: 'black', fontSize: 24, showBackground: true, inFront: true}});
        view.zoomTo();
        view.render();
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
        ballnstick_viz(molecule)
