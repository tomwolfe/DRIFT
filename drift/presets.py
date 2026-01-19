from .topology import Topology
from .metabolic import MetabolicBridge, BridgeBuilder

def get_human_cancer_topology() -> Topology:
    """Returns a signaling topology common in human cancer (includes AMPK)."""
    return Topology(
        species=["PI3K", "AKT", "mTOR", "AMPK"],
        parameters={
            "k_pi3k_base": 0.1,
            "k_pi3k_deg": 0.1,
            "k_akt_act": 0.5,
            "k_akt_deact": 0.1,
            "k_mtor_act": 0.5,
            "k_mtor_deact": 0.1,
            "k_ampk_base": 0.2,
            "k_ampk_deg": 0.1,
        },
        name="Human_Cancer_Signaling"
    )

def get_human_cancer_bridge() -> MetabolicBridge:
    """
    Returns a pre-configured MetabolicBridge for human cancer research.
    Maps PI3K/AKT/mTOR signaling to Recon1 metabolic subsystems.
    """
    builder = BridgeBuilder()
    builder.set_species_names(["PI3K", "AKT", "mTOR", "AMPK"])
    
    # Warburg Effect: mTOR increases glucose uptake and lactate secretion
    builder.add_mapping(protein_name="mTOR", reaction_id="EX_glc__D_e", influence="positive", base_vmax=15.0)
    builder.add_mapping(protein_name="mTOR", reaction_id="EX_lac__L_e", influence="positive", base_vmax=20.0)
    
    # Glutaminolysis: AKT increases glutamine uptake
    builder.add_mapping(protein_name="AKT", reaction_id="EX_gln__L_e", influence="positive", base_vmax=5.0)
    
    # Energy sensing: AMPK activates fatty acid oxidation and oxidative phosphorylation
    builder.add_mapping(protein_name="AMPK", reaction_id="EX_o2_e", influence="positive", base_vmax=20.0)
    
    # Reverse mappings (Metabolism -> Signaling)
    # Low energy (low ATP synthase flux) activates AMPK
    builder.add_reverse_mapping(flux_id="ATPS4r", species_name="AMPK", influence="negative", weight=0.8)
    
    # Growth influences mTOR
    builder.add_reverse_mapping(flux_id="BIOMASS_RECON1", species_name="mTOR", influence="positive", weight=0.5)
    
    return builder.build()
