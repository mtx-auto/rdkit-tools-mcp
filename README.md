# RDKit Tools MCP Server

A Model Context Protocol (MCP) server exposing core [RDKit](https://github.com/rdkit/rdkit) cheminformatics
operations as tools, for use by AI agents (e.g. DIAL toolsets).

Complements other chemistry MCPs (SMILES visualization, ChEMBL database access) by covering the
in-process cheminformatics operations they don't: descriptor/fingerprint calculation, similarity and
substructure search, reaction enumeration, standardization, and 3D conformer generation.

## Tools

| Tool | Description |
|---|---|
| `canonicalize_smiles` | Parse a SMILES string, return its canonical form and basic identity info |
| `standardize_molecule` | Clean up a molecule: strip salts/solvents, keep the largest fragment, neutralize charges |
| `calculate_descriptors` | Compute descriptors (MolWt, LogP, TPSA, HBD/HBA, ring counts, QED, ...) for one SMILES |
| `batch_calculate_descriptors` | Same, for a batch of SMILES |
| `calculate_fingerprint` | Compute a fingerprint (`morgan`, `rdkit`, or `maccs`) and return its on-bit indices |
| `tanimoto_similarity` | Tanimoto similarity between two molecules |
| `similarity_search` | Rank candidate SMILES by similarity to a query molecule |
| `substructure_search` | Check which candidates match a SMARTS substructure pattern |
| `enumerate_reaction` | Apply a reaction SMARTS/SMIRKS template to reactants, return distinct products |
| `generate_3d_conformer` | Generate (optionally force-field optimized) 3D conformers as MOL blocks |

## Running locally

```bash
pip install -r requirements.txt
python server.py --host 127.0.0.1 --port 8080
```

## Docker

```bash
docker build -t rdkit-tools-mcp .
docker run -p 8080:8080 rdkit-tools-mcp
```

The server speaks MCP over HTTP Streamable transport at `/mcp`.
