# ProStructure - 3D

This Streamlit app predicts and visualizes protein structures, providing key analytical insights such as structural quality, residue composition, hydrophobicity, and other metrics.

## Features

* **Interactive 3D Structure Viewer** powered by `py3Dmol` and `stmol`.
* **Quick Statistics Panel** including:

  * Sequence length
  * Molecular weight
  * Estimated isoelectric point (pI)
  * Percentage of residues with high pLDDT confidence scores
  * Solvent Accessible Surface Area (SASA)
  * Disulfide bond count
  * Secondary structure breakdown (helix, sheet, coil)
* **Per-residue Analysis** with hydrophobicity and pLDDT plots.
* **Charge & Composition** visualization with charts.
* **Modern UI** with styled containers and responsive layout.

## Requirements

Install dependencies from `requirements.txt`:

```bash
pip install -r requirements.txt
```

## Usage

Run the app locally:

```bash
streamlit run streamlit_app.py
```

## File Inputs

* The app uses a PDB file (`predicted.pdb`) as input.
* Future updates will support uploading sequences in FASTA format for on-the-fly prediction.

## UI Layout

* **3D Structure**: Left panel with large interactive model.
* **Quick Statistics**: Right panel with structural and sequence data.
* **Charge & Composition**: Visual charts showing residue type distribution and charge.
* **Per-residue Analyses**: Detailed plDDT and hydrophobicity plots.

## Customization

You can adjust:

* Background colors and styles in the CSS section.
* Default coloring scheme for 3D models (`spectrum`, `chain`, `b-factor`, etc.).

## Planned Enhancements

* Automated molecular weight and pI calculation.
* Dynamic pocket and ligand-binding site detection.
* Mutation effect prediction.

---

**Author:** Keeistu M S
**License:** MIT
