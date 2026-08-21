import asyncio
import json
import sys
from server import RDKitToolsMCP


async def main():
    server = RDKitToolsMCP()
    tools = await server.mcp.list_tools()
    names = sorted(t.name for t in tools)
    print("TOOLS:", names)
    assert len(names) == 10, f"expected 10 tools, got {len(names)}"

    async def call(name, **kwargs):
        result = await server.mcp.call_tool(name, kwargs)
        if result.is_error:
            raise RuntimeError(f"tool {name} returned error: {result.content}")
        text = result.content[0].text
        return json.loads(text)

    aspirin = "CC(=O)Oc1ccccc1C(=O)O"
    ibuprofen = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"

    out = await call("canonicalize_smiles", smiles=aspirin)
    print("canonicalize_smiles:", out)
    assert out["valid"]

    out = await call("canonicalize_smiles", smiles="not a smiles")
    print("canonicalize_smiles(invalid):", out)
    assert not out["valid"]

    out = await call("standardize_molecule", smiles="CC(=O)[O-].[Na+]")
    print("standardize_molecule:", out)
    assert "error" not in out

    out = await call("calculate_descriptors", smiles=aspirin)
    print("calculate_descriptors:", out["descriptors"].keys())
    assert "QED" in out["descriptors"]

    out = await call("calculate_descriptors", smiles=aspirin, descriptor_names=["MolWt", "bogus"])
    print("calculate_descriptors(subset):", out)
    assert out["unknown_descriptor_names"] == ["bogus"]

    out = await call("batch_calculate_descriptors", smiles_list=[aspirin, ibuprofen, "bad"])
    print("batch_calculate_descriptors count:", out["count"])
    assert out["count"] == 3

    out = await call("calculate_fingerprint", smiles=aspirin, fp_type="morgan")
    print("calculate_fingerprint(morgan) on_bits:", out["num_on_bits"])
    assert out["num_bits"] == 2048

    out = await call("calculate_fingerprint", smiles=aspirin, fp_type="maccs")
    print("calculate_fingerprint(maccs) num_bits:", out["num_bits"])

    out = await call("tanimoto_similarity", smiles_a=aspirin, smiles_b=ibuprofen)
    print("tanimoto_similarity:", out["tanimoto_similarity"])
    assert 0.0 <= out["tanimoto_similarity"] <= 1.0

    out = await call("similarity_search", query_smiles=aspirin, candidate_smiles=[aspirin, ibuprofen, "bad"], top_k=5)
    print("similarity_search:", out["results"])
    assert out["results"][0]["smiles"] == aspirin
    assert out["errors"][0]["smiles"] == "bad"

    out = await call("substructure_search", smarts="c1ccccc1", candidate_smiles=[aspirin, "CCO"])
    print("substructure_search:", out["results"])
    assert out["results"][0]["matched"] is True
    assert out["results"][1]["matched"] is False

    # Esterification: acid + alcohol -> ester
    rxn_smarts = "[C:1](=[O:2])[OH:3].[OH:4][C:5]>>[C:1](=[O:2])[O:4][C:5]"
    out = await call("enumerate_reaction", reaction_smarts=rxn_smarts, reactant_smiles=["CC(=O)O", "CCO"])
    print("enumerate_reaction:", out)
    assert out["products"], "expected at least one product"

    out = await call("generate_3d_conformer", smiles="CCO", num_conformers=2, optimize=True)
    print("generate_3d_conformer:", out["num_conformers_generated"], "energies:", [c["energy"] for c in out["conformers"]])
    assert out["num_conformers_generated"] == 2
    assert all(c["mol_block"] for c in out["conformers"])

    print("\nALL CHECKS PASSED")


if __name__ == "__main__":
    asyncio.run(main())
