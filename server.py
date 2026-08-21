#!/usr/bin/env python3
"""
RDKit Tools MCP Server
A Model Context Protocol server exposing core RDKit cheminformatics tools:
canonicalization/standardization, descriptors, fingerprints, similarity
search, substructure search, reaction enumeration, and 3D conformer
generation.
"""

import asyncio
import argparse
import json
import logging
import os
import sys
from typing import Any, Dict, List, Optional

from mcp.server.mcpserver import MCPServer
from mcp.server.transport_security import TransportSecuritySettings

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

try:
    from rdkit import Chem, RDLogger
    from rdkit.Chem import AllChem, Descriptors, Crippen, rdMolDescriptors, MACCSkeys, QED
    from rdkit.Chem.MolStandardize import rdMolStandardize
    from rdkit import DataStructs
    RDKIT_AVAILABLE = True
    RDLogger.DisableLog("rdApp.*")
except ImportError as e:
    RDKIT_AVAILABLE = False
    logger.warning(f"RDKit not available: {e}")

DESCRIPTOR_FUNCS = {
    "MolWt": Descriptors.MolWt,
    "ExactMolWt": rdMolDescriptors.CalcExactMolWt,
    "MolLogP": Crippen.MolLogP,
    "TPSA": rdMolDescriptors.CalcTPSA,
    "NumHDonors": rdMolDescriptors.CalcNumHBD,
    "NumHAcceptors": rdMolDescriptors.CalcNumHBA,
    "NumRotatableBonds": rdMolDescriptors.CalcNumRotatableBonds,
    "RingCount": rdMolDescriptors.CalcNumRings,
    "NumAromaticRings": rdMolDescriptors.CalcNumAromaticRings,
    "FractionCSP3": rdMolDescriptors.CalcFractionCSP3,
    "HeavyAtomCount": lambda mol: mol.GetNumHeavyAtoms(),
    "NumAtoms": lambda mol: mol.GetNumAtoms(),
    "MolecularFormula": rdMolDescriptors.CalcMolFormula,
    "NumRadicalElectrons": Descriptors.NumRadicalElectrons,
    "NumValenceElectrons": Descriptors.NumValenceElectrons,
    "QED": QED.qed,
} if RDKIT_AVAILABLE else {}


def _require_rdkit():
    if not RDKIT_AVAILABLE:
        raise RuntimeError("RDKit is not available in this environment")


def _mol_from_smiles(smiles: str) -> "Chem.Mol":
    _require_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")
    return mol


def _get_fingerprint(mol: "Chem.Mol", fp_type: str, radius: int, n_bits: int):
    fp_type = fp_type.lower()
    if fp_type == "morgan":
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)
    if fp_type == "rdkit":
        return Chem.RDKFingerprint(mol, fpSize=n_bits)
    if fp_type == "maccs":
        return MACCSkeys.GenMACCSKeys(mol)
    raise ValueError(f"Unknown fp_type '{fp_type}'. Use 'morgan', 'rdkit', or 'maccs'.")


class RDKitToolsMCP:
    """RDKit cheminformatics MCP server implementation using mcp.server.mcpserver.MCPServer"""

    def __init__(self):
        self.mcp = MCPServer("RDKit Tools")
        self.setup_tools()

    def setup_tools(self):
        @self.mcp.tool()
        async def canonicalize_smiles(smiles: str) -> str:
            """Parse a SMILES string and return its canonical form plus basic identity info"""
            try:
                mol = _mol_from_smiles(smiles)
                return json.dumps({
                    "valid": True,
                    "input_smiles": smiles,
                    "canonical_smiles": Chem.MolToSmiles(mol, canonical=True),
                    "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
                    "num_atoms": mol.GetNumAtoms(),
                    "num_heavy_atoms": mol.GetNumHeavyAtoms(),
                })
            except Exception as e:
                return json.dumps({"valid": False, "error": str(e)})

        @self.mcp.tool()
        async def standardize_molecule(smiles: str) -> str:
            """Standardize a molecule: clean up, strip salts/solvents (keep largest fragment), and neutralize charges"""
            try:
                mol = _mol_from_smiles(smiles)
                cleaned = rdMolStandardize.Cleanup(mol)
                parent = rdMolStandardize.FragmentParent(cleaned)
                uncharger = rdMolStandardize.Uncharger()
                neutral = uncharger.uncharge(parent)
                return json.dumps({
                    "input_smiles": smiles,
                    "standardized_smiles": Chem.MolToSmiles(neutral, canonical=True),
                })
            except Exception as e:
                return json.dumps({"error": str(e)})

        @self.mcp.tool()
        async def calculate_descriptors(smiles: str, descriptor_names: Optional[List[str]] = None) -> str:
            """Calculate molecular descriptors (MolWt, LogP, TPSA, HBD/HBA, RingCount, QED, etc.) for a SMILES string.
            Pass descriptor_names to restrict to a subset; omit for the full default set.
            """
            try:
                mol = _mol_from_smiles(smiles)
                names = descriptor_names or list(DESCRIPTOR_FUNCS.keys())
                result = {}
                unknown = []
                for name in names:
                    func = DESCRIPTOR_FUNCS.get(name)
                    if func is None:
                        unknown.append(name)
                        continue
                    result[name] = func(mol)
                response = {"smiles": smiles, "descriptors": result}
                if unknown:
                    response["unknown_descriptor_names"] = unknown
                    response["available_descriptor_names"] = list(DESCRIPTOR_FUNCS.keys())
                return json.dumps(response, indent=2)
            except Exception as e:
                return json.dumps({"error": str(e)})

        @self.mcp.tool()
        async def batch_calculate_descriptors(smiles_list: List[str], descriptor_names: Optional[List[str]] = None) -> str:
            """Calculate molecular descriptors for a batch of SMILES strings, one result entry per input"""
            results = []
            for smiles in smiles_list:
                parsed = json.loads(await calculate_descriptors(smiles, descriptor_names))
                parsed["smiles"] = smiles
                results.append(parsed)
            return json.dumps({"count": len(results), "results": results}, indent=2)

        @self.mcp.tool()
        async def calculate_fingerprint(smiles: str, fp_type: str = "morgan", radius: int = 2, n_bits: int = 2048) -> str:
            """Calculate a molecular fingerprint ('morgan', 'rdkit', or 'maccs') and return its on-bit indices"""
            try:
                mol = _mol_from_smiles(smiles)
                fp = _get_fingerprint(mol, fp_type, radius, n_bits)
                on_bits = list(fp.GetOnBits())
                return json.dumps({
                    "smiles": smiles,
                    "fp_type": fp_type,
                    "num_bits": fp.GetNumBits(),
                    "num_on_bits": len(on_bits),
                    "on_bits": on_bits,
                })
            except Exception as e:
                return json.dumps({"error": str(e)})

        @self.mcp.tool()
        async def tanimoto_similarity(smiles_a: str, smiles_b: str, fp_type: str = "morgan", radius: int = 2, n_bits: int = 2048) -> str:
            """Compute the Tanimoto similarity between two molecules given as SMILES"""
            try:
                mol_a = _mol_from_smiles(smiles_a)
                mol_b = _mol_from_smiles(smiles_b)
                fp_a = _get_fingerprint(mol_a, fp_type, radius, n_bits)
                fp_b = _get_fingerprint(mol_b, fp_type, radius, n_bits)
                similarity = DataStructs.TanimotoSimilarity(fp_a, fp_b)
                return json.dumps({
                    "smiles_a": smiles_a,
                    "smiles_b": smiles_b,
                    "fp_type": fp_type,
                    "tanimoto_similarity": similarity,
                })
            except Exception as e:
                return json.dumps({"error": str(e)})

        @self.mcp.tool()
        async def similarity_search(query_smiles: str, candidate_smiles: List[str], fp_type: str = "morgan",
                                     radius: int = 2, n_bits: int = 2048, top_k: int = 10) -> str:
            """Rank a list of candidate SMILES by Tanimoto similarity to a query molecule, most similar first"""
            try:
                query_mol = _mol_from_smiles(query_smiles)
                query_fp = _get_fingerprint(query_mol, fp_type, radius, n_bits)
                scored = []
                errors = []
                for candidate in candidate_smiles:
                    try:
                        mol = _mol_from_smiles(candidate)
                        fp = _get_fingerprint(mol, fp_type, radius, n_bits)
                        scored.append({
                            "smiles": candidate,
                            "similarity": DataStructs.TanimotoSimilarity(query_fp, fp),
                        })
                    except Exception as e:
                        errors.append({"smiles": candidate, "error": str(e)})
                scored.sort(key=lambda x: x["similarity"], reverse=True)
                response = {"query_smiles": query_smiles, "fp_type": fp_type, "results": scored[:top_k]}
                if errors:
                    response["errors"] = errors
                return json.dumps(response, indent=2)
            except Exception as e:
                return json.dumps({"error": str(e)})

        @self.mcp.tool()
        async def substructure_search(smarts: str, candidate_smiles: List[str]) -> str:
            """Check which candidate SMILES contain a substructure matching the given SMARTS pattern"""
            try:
                _require_rdkit()
                pattern = Chem.MolFromSmarts(smarts)
                if pattern is None:
                    raise ValueError(f"Invalid SMARTS pattern: {smarts}")
                matches = []
                errors = []
                for candidate in candidate_smiles:
                    try:
                        mol = _mol_from_smiles(candidate)
                        match_indices = mol.GetSubstructMatches(pattern)
                        matches.append({
                            "smiles": candidate,
                            "matched": len(match_indices) > 0,
                            "num_matches": len(match_indices),
                            "match_atom_indices": [list(m) for m in match_indices],
                        })
                    except Exception as e:
                        errors.append({"smiles": candidate, "error": str(e)})
                response = {"smarts": smarts, "results": matches}
                if errors:
                    response["errors"] = errors
                return json.dumps(response, indent=2)
            except Exception as e:
                return json.dumps({"error": str(e)})

        @self.mcp.tool()
        async def enumerate_reaction(reaction_smarts: str, reactant_smiles: List[str]) -> str:
            """Apply a reaction SMARTS/SMIRKS template to a set of reactant SMILES and return the distinct product SMILES"""
            try:
                _require_rdkit()
                rxn = AllChem.ReactionFromSmarts(reaction_smarts)
                if rxn is None:
                    raise ValueError(f"Invalid reaction SMARTS: {reaction_smarts}")
                expected = rxn.GetNumReactantTemplates()
                if len(reactant_smiles) != expected:
                    raise ValueError(
                        f"Reaction expects {expected} reactant(s), got {len(reactant_smiles)}"
                    )
                reactant_mols = tuple(_mol_from_smiles(s) for s in reactant_smiles)
                product_sets = rxn.RunReactants(reactant_mols)
                products = set()
                for product_set in product_sets:
                    for product in product_set:
                        try:
                            Chem.SanitizeMol(product)
                            products.add(Chem.MolToSmiles(product, canonical=True))
                        except Exception:
                            continue
                return json.dumps({
                    "reaction_smarts": reaction_smarts,
                    "reactant_smiles": reactant_smiles,
                    "num_product_sets": len(product_sets),
                    "products": sorted(products),
                })
            except Exception as e:
                return json.dumps({"error": str(e)})

        @self.mcp.tool()
        async def generate_3d_conformer(smiles: str, num_conformers: int = 1, optimize: bool = True, seed: int = 42) -> str:
            """Generate 3D conformer(s) for a SMILES string, optionally force-field optimized (MMFF, falling back to UFF), returned as MOL blocks"""
            try:
                mol = _mol_from_smiles(smiles)
                mol = Chem.AddHs(mol)
                conf_ids = list(AllChem.EmbedMultipleConfs(
                    mol, numConfs=max(1, num_conformers), randomSeed=seed
                ))
                if not conf_ids:
                    raise ValueError("Embedding failed to generate any conformers")

                energies = {cid: None for cid in conf_ids}
                if optimize:
                    try:
                        mmff_results = AllChem.MMFFOptimizeMoleculeConfs(mol)
                        for cid, (_, energy) in zip(conf_ids, mmff_results):
                            energies[cid] = energy
                    except Exception:
                        try:
                            uff_results = AllChem.UFFOptimizeMoleculeConfs(mol)
                            for cid, (_, energy) in zip(conf_ids, uff_results):
                                energies[cid] = energy
                        except Exception:
                            pass

                conformers = [
                    {
                        "conformer_id": cid,
                        "energy": energies.get(cid),
                        "mol_block": Chem.MolToMolBlock(mol, confId=cid),
                    }
                    for cid in conf_ids
                ]
                return json.dumps({
                    "smiles": smiles,
                    "num_conformers_requested": num_conformers,
                    "num_conformers_generated": len(conformers),
                    "optimized": optimize,
                    "conformers": conformers,
                })
            except Exception as e:
                return json.dumps({"error": str(e)})


async def run_http_streamable_server(server: RDKitToolsMCP, host: str, port: int):
    """Run the MCP server using HTTP Streamable transport"""
    logger.info(f"Starting RDKit Tools MCP Server with HTTP Streamable transport on {host}:{port}")

    import uvicorn
    app = server.mcp.streamable_http_app(
        host=host,
        transport_security=TransportSecuritySettings(enable_dns_rebinding_protection=False),
    )
    config = uvicorn.Config(app, host=host, port=port, log_level="info")
    uv_server = uvicorn.Server(config)
    await uv_server.serve()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="RDKit Tools MCP Server - cheminformatics tools over the Model Context Protocol",
    )
    parser.add_argument("--host", "-H", default=os.getenv("MCP_HOST", "127.0.0.1"))
    parser.add_argument("--port", "-p", type=int, default=int(os.getenv("MCP_PORT", "8080")))
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--version", action="version", version="RDKit Tools MCP Server 1.0.0")
    return parser.parse_args()


async def main():
    args = parse_arguments()
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    logger.info("RDKit Tools MCP Server Configuration:")
    logger.info(f"  Host: {args.host}")
    logger.info(f"  Port: {args.port}")
    logger.info(f"  RDKit Available: {RDKIT_AVAILABLE}")

    try:
        server = RDKitToolsMCP()
        await run_http_streamable_server(server, args.host, args.port)
    except KeyboardInterrupt:
        logger.info("Server stopped by user")
    except Exception as e:
        logger.error(f"Server error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    asyncio.run(main())
