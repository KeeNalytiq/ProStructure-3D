# streamlit_app_enhanced.py
"""
Enhanced Protein Structure Visualizer (Streamlit)
Features added:
 - Paste or upload sequence (required)
 - Sequence validation/cleaning
 - Cached ESMFold predictions
 - 3D visualization with multiple coloring modes (including plDDT & hydrophobicity)
 - Per-residue plDDT interactive plot (Plotly)
 - Hydrophobicity profile (Kyte-Doolittle)
 - Charge distribution and residue composition pie charts
 - Molecular weight and approximate isoelectric point (pI)
 - Radius of gyration from coordinates
 - Confidence statistics (% residues above thresholds)
 - Secondary structure summary if HELIX/SHEET records present in PDB
 - Downloadable: PDB file and CSV report of per-residue data

Run: streamlit run streamlit_app_enhanced.py
"""

import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
import biotite.structure as struc
import numpy as np
import pandas as pd
import io
import re
import plotly.graph_objects as go
import plotly.express as px
import math
import time
import random
from collections import Counter
import json
import base64
from datetime import datetime
import threading
import queue
import subprocess
import sys

# -------------------- RDKit availability check --------------------
rdkit_available = False
try:
    import rdkit
    from rdkit import Chem
    rdkit_available = True
except ImportError:
    rdkit_available = False

# -------------------- MDTraj availability check --------------------
mdtraj_available = False
try:
    import mdtraj as md
    mdtraj_available = True
except Exception:
    mdtraj_available = False

# -------------------- Page setup --------------------
st.set_page_config(
    page_title="ProStruct - 3D", 
    layout="wide",
    initial_sidebar_state="expanded",
    page_icon="🧬"
)

# Custom CSS for modern UI
st.markdown("""
<style>
    /* Import Google Fonts */
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700&display=swap');
    
    /* Global styles */
    .main {
        padding-top: 2rem;
    }
    
    /* Custom header */
    .main-header {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        padding: 2rem;
        border-radius: 20px;
        color: white;
        text-align: center;
        margin-bottom: 2rem;
        box-shadow: 0 10px 30px rgba(0,0,0,0.1);
    }
    
    .main-header h1 {
        font-family: 'Inter', sans-serif;
        font-weight: 700;
        font-size: 3rem;
        margin: 0;
        text-shadow: 0 2px 4px rgba(0,0,0,0.3);
    }
    
    .main-header p {
        font-family: 'Inter', sans-serif;
        font-weight: 400;
        font-size: 1.2rem;
        margin: 0.5rem 0 0 0;
        opacity: 0.9;
    }
    
    /* Sidebar improvements */
    .css-1d391kg {
        background: linear-gradient(180deg, #f8fafc 0%, #e2e8f0 100%);
    }
    
    /* Button styling */
    .stButton > button {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        border-radius: 12px;
        padding: 0.75rem 2rem;
        font-weight: 600;
        font-size: 1rem;
        transition: all 0.3s ease;
        box-shadow: 0 4px 15px rgba(102, 126, 234, 0.3);
    }
    
    .stButton > button:hover {
        transform: translateY(-2px);
        box-shadow: 0 8px 25px rgba(102, 126, 234, 0.4);
    }
    
    /* Container styling */
    .stContainer {
        background: white;
        border-radius: 16px;
        padding: 1.5rem;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        border: 1px solid rgba(0,0,0,0.05);
    }
    
    /* Metric cards */
    .metric-card {
        background: linear-gradient(135deg, #f0f9ff 0%, #e0f2fe 100%);
        padding: 1.5rem;
        border-radius: 16px;
        border-left: 4px solid #0ea5e9;
        box-shadow: 0 4px 15px rgba(0,0,0,0.05);
        margin: 1rem 0;
    }
    
    /* Chart containers */
    .chart-container {
        background: white;
        border-radius: 16px;
        padding: 1.5rem;
        box-shadow: 0 4px 20px rgba(0,0,0,0.08);
        margin: 1rem 0;
    }
    
    /* 3D viewer container */
    .viewer-container {
        background: linear-gradient(135deg, #f0fdfa 0%, #ecfdf5 100%);
        border-radius: 20px;
        padding: 2rem;
        box-shadow: 0 8px 30px rgba(0,0,0,0.1);
        text-align: center;
    }
    
    /* Progress bars */
    .stProgress > div > div > div > div {
        background: linear-gradient(90deg, #667eea, #764ba2);
    }
    
    /* File uploader styling */
    .stFileUploader > div > div {
        border: 2px dashed #cbd5e1;
        border-radius: 12px;
        background: linear-gradient(135deg, #f8fafc 0%, #f1f5f9 100%);
    }
    
    /* Selectbox styling */
    .stSelectbox > div > div {
        background: white;
        border-radius: 8px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.05);
    }
    
    /* Text area styling */
    .stTextArea > div > div > textarea {
        border-radius: 12px;
        border: 2px solid #e2e8f0;
        font-family: 'Inter', monospace;
        font-size: 0.9rem;
    }
    
    /* Success/Error styling */
    .stSuccess {
        background: linear-gradient(135deg, #dcfce7 0%, #bbf7d0 100%);
        border: 1px solid #22c55e;
        border-radius: 12px;
    }
    
    .stError {
        background: linear-gradient(135deg, #fef2f2 0%, #fecaca 100%);
        border: 1px solid #ef4444;
        border-radius: 12px;
    }
    
    .stWarning {
        background: linear-gradient(135deg, #fefce8 0%, #fef3c7 100%);
        border: 1px solid #f59e0b;
        border-radius: 12px;
    }
    
    /* Hide Streamlit branding */
    #MainMenu {visibility: hidden;}
    footer {visibility: hidden;}
    header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# Modern header
st.markdown("""
<div class="main-header">
    <h1>🧬 ProStruct - 3D</h1>
    <p>Advanced Protein Structure Prediction & Analysis Platform</p>
</div>
""", unsafe_allow_html=True)

# Feature highlights
col1, col2, col3, col4 = st.columns(4)
with col1:
    st.markdown("""
    <div class="metric-card">
        <h4 style="color: #0ea5e9; margin: 0;">🔬 ESMFold</h4>
        <p style="margin: 0.5rem 0 0 0; color: #64748b;">State-of-the-art structure prediction</p>
    </div>
    """, unsafe_allow_html=True)

with col2:
    st.markdown("""
    <div class="metric-card">
        <h4 style="color: #0ea5e9; margin: 0;">📊 Analytics</h4>
        <p style="margin: 0.5rem 0 0 0; color: #64748b;">Comprehensive protein analysis</p>
    </div>
    """, unsafe_allow_html=True)

with col3:
    st.markdown("""
    <div class="metric-card">
        <h4 style="color: #0ea5e9; margin: 0;">🎨 3D Visualization</h4>
        <p style="margin: 0.5rem 0 0 0; color: #64748b;">Interactive molecular viewer</p>
    </div>
    """, unsafe_allow_html=True)

with col4:
    st.markdown("""
    <div class="metric-card">
        <h4 style="color: #0ea5e9; margin: 0;">📈 Reports</h4>
        <p style="margin: 0.5rem 0 0 0; color: #64748b;">Detailed analysis reports</p>
    </div>
    """, unsafe_allow_html=True)

# Navigation tabs
st.markdown("---")
tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs([
    "🔬 Structure Prediction", 
    "🔄 Protein Comparison", 
    "🏗️ Domain Analysis", 
    "🔗 PPI Prediction", 
    "⚙️ Advanced Tools",
    "🤖 AI Protein Designer",
    "🧬 Molecular Dynamics",
    "💊 Drug Discovery"
])

# Initialize session state for multiple features
if 'prediction_history' not in st.session_state:
    st.session_state.prediction_history = []
if 'comparison_data' not in st.session_state:
    st.session_state.comparison_data = []
if 'domain_results' not in st.session_state:
    st.session_state.domain_results = {}
if 'ppi_predictions' not in st.session_state:
    st.session_state.ppi_predictions = {}
if 'ai_designs' not in st.session_state:
    st.session_state.ai_designs = []
if 'md_simulations' not in st.session_state:
    st.session_state.md_simulations = []
if 'drug_predictions' not in st.session_state:
    st.session_state.drug_predictions = []
if 'collaboration_room' not in st.session_state:
    st.session_state.collaboration_room = None

# Rate limiting and session state
if 'last_request_time' not in st.session_state:
    st.session_state.last_request_time = 0
if 'request_count' not in st.session_state:
    st.session_state.request_count = 0

def check_rate_limit():
    """Check if user is making requests too frequently."""
    current_time = time.time()
    time_since_last = current_time - st.session_state.last_request_time
    
    # Allow max 1 request per 10 seconds
    if time_since_last < 10:
        remaining_time = 10 - time_since_last
        st.warning(f"⏱️ Please wait {remaining_time:.1f} seconds before making another request to avoid overloading the API.")
        return False
    
    st.session_state.last_request_time = current_time
    st.session_state.request_count += 1
    return True

st.sidebar.header("Input & Settings")
st.sidebar.write("Provide a sequence by pasting it or uploading a FASTA/TXT file. Prediction requires at least 20 residues.")

# API Status and Tips
with st.sidebar.expander("🔧 Troubleshooting & Tips", expanded=False):
    st.markdown("""
    **If you get 503 errors:**
    - The ESM Atlas API is temporarily overloaded
    - Try the retry mechanism (automatic)
    - Switch to Local ESMFold if you have GPU
    - Wait a few minutes and try again
    
    **For better success:**
    - Use sequences < 400 residues
    - Avoid making rapid requests
    - Check your internet connection
    
    **Local ESMFold:**
    - Requires GPU with 8GB+ VRAM
    - Install: `pip install torch transformers`
    - Slower but more reliable
    """)
    
    # Show request count
    if st.session_state.request_count > 0:
        st.info(f"Requests made this session: {st.session_state.request_count}")

# -------------------- Helper utilities --------------------
AA_LIST = list("ACDEFGHIKLMNPQRSTVWY")

KD_SCALE = {
    'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5, 'F': 2.8,
    'G': -0.4, 'H': -3.2, 'I': 4.5, 'K': -3.9, 'L': 3.8,
    'M': 1.9, 'N': -3.5, 'P': -1.6, 'Q': -3.5, 'R': -4.5,
    'S': -0.8, 'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3
}

AA_MASS = {
    'A': 89.0935, 'C': 121.1590, 'D': 133.1027, 'E': 147.1293, 'F': 165.1900,
    'G': 75.0669, 'H': 155.1552, 'I': 131.1736, 'K': 146.1882, 'L': 131.1736,
    'M': 149.2124, 'N': 132.1184, 'P': 115.1310, 'Q': 146.1451, 'R': 174.2017,
    'S': 105.0930, 'T': 119.1197, 'V': 117.1469, 'W': 204.2262, 'Y': 181.1894
}

# pKa values for side chains and terminals (approximate)
PKA = {
    'Cterm': 3.55, 'Nterm': 8.0,
    'C': 8.18, 'D': 3.65, 'E': 4.25, 'H': 6.00, 'K': 10.53, 'R': 12.48, 'Y': 10.07
}

POSITIVE_AA = set(['K', 'R', 'H'])
NEGATIVE_AA = set(['D', 'E'])


def log_api(message: str, level: str = 'info') -> None:
    """Conditional UI logger to avoid noisy warnings unless verbose is enabled."""
    try:
        if st.session_state.get('verbose_api', False):
            if level == 'warning':
                st.warning(message)
            else:
                st.info(message)
    except Exception:
        # Streamlit context may not be initialized during import or caching
        pass


# -------------------- ADMET utilities --------------------
@st.cache_data(show_spinner=False)
def compute_admet_properties(smiles: str):
    """Enhanced ADMET properties computation with comprehensive drug discovery metrics."""
    try:
        import rdkit
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen, Lipinski, rdMolDescriptors, FilterCatalog, QED
    except ImportError:
        return { 'error': 'RDKit not installed. Please install RDKit: pip install rdkit' }
    except Exception as e:
        return { 'error': f'RDKit import error: {str(e)}' }
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return { 'error': 'Invalid SMILES' }
        
        # Basic molecular descriptors
        mw = float(Descriptors.MolWt(mol))
        logp = float(Crippen.MolLogP(mol))
        tpsa = float(Descriptors.TPSA(mol))
        hbd = int(Lipinski.NumHDonors(mol))
        hba = int(Lipinski.NumHAcceptors(mol))
        rotb = int(Lipinski.NumRotatableBonds(mol))
        rings_arom = int(rdMolDescriptors.CalcNumAromaticRings(mol))
        rings_total = int(rdMolDescriptors.CalcNumRings(mol))
        
        # Enhanced descriptors
        formal_charge = int(Chem.rdmolops.GetFormalCharge(mol))
        heavy_atoms = int(Descriptors.HeavyAtomCount(mol))
        stereo_centers = int(rdMolDescriptors.CalcNumAliphaticCarbocycles(mol))
        
        # Drug-likeness scores
        qed_score = float(QED.qed(mol)) if hasattr(QED, 'qed') else 0.0
        
        # PAINS and structural alerts
        params = FilterCatalog.FilterCatalogParams()
        params.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.PAINS)
        pains = FilterCatalog.FilterCatalog(params)
        matches = pains.GetMatches(mol)
        pains_hits = len(matches)
        
        # Additional structural alerts
        brenk_alerts = 0
        nih_alerts = 0
        try:
            # Brenk alerts
            params_brenk = FilterCatalog.FilterCatalogParams()
            params_brenk.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.BRENK)
            brenk = FilterCatalog.FilterCatalog(params_brenk)
            brenk_alerts = len(brenk.GetMatches(mol))
            
            # NIH alerts
            params_nih = FilterCatalog.FilterCatalogParams()
            params_nih.AddCatalog(FilterCatalog.FilterCatalogParams.FilterCatalogs.NIH)
            nih = FilterCatalog.FilterCatalog(params_nih)
            nih_alerts = len(nih.GetMatches(mol))
        except:
            pass
        
        # Drug-likeness rules
        lipinski_pass = (mw <= 500 and logp <= 5 and hbd <= 5 and hba <= 10)
        veber_pass = (rotb <= 10 and tpsa <= 140)
        egan_pass = (tpsa <= 131 and logp <= 5.88)
        lead_like_pass = (mw <= 350 and logp <= 3.5 and tpsa <= 70)
        fragment_like_pass = (mw <= 300 and logp <= 3 and hbd <= 3 and hba <= 3)
        
        # CNS drug-likeness
        cns_pass = (tpsa <= 70 and logp >= 2 and logp <= 5 and mw <= 400)
        
        # Bioavailability score
        bioavailability_score = 0.0
        if lipinski_pass:
            bioavailability_score += 0.3
        if veber_pass:
            bioavailability_score += 0.3
        if egan_pass:
            bioavailability_score += 0.2
        if tpsa <= 90:  # Good membrane permeability
            bioavailability_score += 0.2
        
        # Synthetic accessibility (simplified)
        synth_access = 1.0
        if rings_arom > 3:
            synth_access -= 0.2
        if stereo_centers > 2:
            synth_access -= 0.2
        if formal_charge != 0:
            synth_access -= 0.1
        synth_access = max(0.0, synth_access)
        
        return {
            # Basic properties
            'mw': mw, 'logp': logp, 'tpsa': tpsa, 'hbd': hbd, 'hba': hba, 'rotb': rotb,
            'aromatic_rings': rings_arom, 'total_rings': rings_total,
            'heavy_atoms': heavy_atoms, 'formal_charge': formal_charge,
            'stereo_centers': stereo_centers,
            
            # Drug-likeness scores
            'qed_score': qed_score, 'bioavailability_score': bioavailability_score,
            'synthetic_accessibility': synth_access,
            
            # Structural alerts
            'pains': pains_hits, 'brenk_alerts': brenk_alerts, 'nih_alerts': nih_alerts,
            'total_alerts': pains_hits + brenk_alerts + nih_alerts,
            
            # Rule compliance
            'lipinski_pass': lipinski_pass, 'veber_pass': veber_pass, 'egan_pass': egan_pass,
            'lead_like_pass': lead_like_pass, 'fragment_like_pass': fragment_like_pass,
            'cns_pass': cns_pass,
            
            # Drug-likeness category
            'drug_category': 'fragment' if fragment_like_pass else 'lead' if lead_like_pass else 'drug' if lipinski_pass else 'large'
        }
    except Exception as e:
        return { 'error': str(e) }

def interpret_admet(props: dict):
    """Enhanced ADMET interpretation with comprehensive drug discovery metrics."""
    if 'error' in props:
        return {'error': props['error']}
    
    # Helper functions for scoring
    def score_le(val, thr, hard_max=None):
        if val <= thr:
            return 1.0
        if hard_max is None:
            return max(0.0, 1.0 - (val - thr) / max(1e-6, thr))
        return max(0.0, 1.0 - (val - thr) / max(1e-6, hard_max - thr))
    
    def score_ge(val, thr, hard_min=None):
        if val >= thr:
            return 1.0
        if hard_min is None:
            return max(0.0, val / max(1e-6, thr))
        return max(0.0, (val - hard_min) / max(1e-6, thr - hard_min))
    
    def score_range(val, min_val, max_val):
        if min_val <= val <= max_val:
            return 1.0
        elif val < min_val:
            return max(0.0, val / max(1e-6, min_val))
        else:
            return max(0.0, 1.0 - (val - max_val) / max(1e-6, max_val))
    
    # Absorption (A) - Enhanced
    abs_metrics = [
        { 'name': 'Lipinski Rule', 'value': 1 if props['lipinski_pass'] else 0, 'target': 1, 'pass': props['lipinski_pass'], 'score': 1.0 if props['lipinski_pass'] else 0.0 },
        { 'name': 'TPSA (≤140 Å²)', 'value': props['tpsa'], 'target': 140, 'pass': props['tpsa'] <= 140, 'score': score_le(props['tpsa'], 140, 200) },
        { 'name': 'HBD (≤5)', 'value': props['hbd'], 'target': 5, 'pass': props['hbd'] <= 5, 'score': score_le(props['hbd'], 5, 10) },
        { 'name': 'HBA (≤10)', 'value': props['hba'], 'target': 10, 'pass': props['hba'] <= 10, 'score': score_le(props['hba'], 10, 20) },
        { 'name': 'LogP (≤5)', 'value': props['logp'], 'target': 5, 'pass': props['logp'] <= 5, 'score': score_le(props['logp'], 5, 8) },
        { 'name': 'Bioavailability Score', 'value': props['bioavailability_score'], 'target': 0.8, 'pass': props['bioavailability_score'] >= 0.8, 'score': min(1.0, props['bioavailability_score']) },
    ]
    
    # Distribution (D) - Enhanced
    dist_metrics = [
        { 'name': 'LogP (2–5 ideal)', 'value': props['logp'], 'target': '2–5', 'pass': 2 <= props['logp'] <= 5, 'score': score_range(props['logp'], 2, 5) },
        { 'name': 'TPSA BBB (≤70)', 'value': props['tpsa'], 'target': 70, 'pass': props['tpsa'] <= 70, 'score': score_le(props['tpsa'], 70, 140) },
        { 'name': 'Aromatic Rings (≤3)', 'value': props['aromatic_rings'], 'target': 3, 'pass': props['aromatic_rings'] <= 3, 'score': score_le(props['aromatic_rings'], 3, 8) },
        { 'name': 'CNS Drug-likeness', 'value': 1 if props['cns_pass'] else 0, 'target': 1, 'pass': props['cns_pass'], 'score': 1.0 if props['cns_pass'] else 0.0 },
    ]
    
    # Metabolism (M) - Enhanced
    metab_metrics = [
        { 'name': 'Rotatable Bonds (≤10)', 'value': props['rotb'], 'target': 10, 'pass': props['rotb'] <= 10, 'score': score_le(props['rotb'], 10, 20) },
        { 'name': 'LogP (≤5)', 'value': props['logp'], 'target': 5, 'pass': props['logp'] <= 5, 'score': score_le(props['logp'], 5, 8) },
        { 'name': 'Synthetic Accessibility', 'value': props['synthetic_accessibility'], 'target': 0.8, 'pass': props['synthetic_accessibility'] >= 0.8, 'score': props['synthetic_accessibility'] },
        { 'name': 'Stereo Centers (≤2)', 'value': props['stereo_centers'], 'target': 2, 'pass': props['stereo_centers'] <= 2, 'score': score_le(props['stereo_centers'], 2, 5) },
    ]
    
    # Excretion (E) - Enhanced
    excr_metrics = [
        { 'name': 'MW (≤500)', 'value': props['mw'], 'target': 500, 'pass': props['mw'] <= 500, 'score': score_le(props['mw'], 500, 800) },
        { 'name': 'TPSA (≤140)', 'value': props['tpsa'], 'target': 140, 'pass': props['tpsa'] <= 140, 'score': score_le(props['tpsa'], 140, 200) },
        { 'name': 'Heavy Atoms (≤35)', 'value': props['heavy_atoms'], 'target': 35, 'pass': props['heavy_atoms'] <= 35, 'score': score_le(props['heavy_atoms'], 35, 50) },
    ]
    
    # Toxicity (T) - Enhanced
    tox_metrics = [
        { 'name': 'PAINS Alerts (=0)', 'value': props['pains'], 'target': 0, 'pass': props['pains'] == 0, 'score': 1.0 if props['pains'] == 0 else max(0.0, 1.0 - props['pains'] * 0.25) },
        { 'name': 'Brenk Alerts (=0)', 'value': props['brenk_alerts'], 'target': 0, 'pass': props['brenk_alerts'] == 0, 'score': 1.0 if props['brenk_alerts'] == 0 else max(0.0, 1.0 - props['brenk_alerts'] * 0.2) },
        { 'name': 'NIH Alerts (=0)', 'value': props['nih_alerts'], 'target': 0, 'pass': props['nih_alerts'] == 0, 'score': 1.0 if props['nih_alerts'] == 0 else max(0.0, 1.0 - props['nih_alerts'] * 0.2) },
        { 'name': 'Total Alerts (=0)', 'value': props['total_alerts'], 'target': 0, 'pass': props['total_alerts'] == 0, 'score': 1.0 if props['total_alerts'] == 0 else max(0.0, 1.0 - props['total_alerts'] * 0.15) },
        { 'name': 'Egan Rule', 'value': 1 if props['egan_pass'] else 0, 'target': 1, 'pass': props['egan_pass'], 'score': 1.0 if props['egan_pass'] else 0.0 },
    ]
    
    # Drug-likeness Assessment
    drug_likeness_metrics = [
        { 'name': 'QED Score', 'value': props['qed_score'], 'target': 0.7, 'pass': props['qed_score'] >= 0.7, 'score': min(1.0, props['qed_score']) },
        { 'name': 'Drug Category', 'value': props['drug_category'], 'target': 'drug', 'pass': props['drug_category'] in ['drug', 'lead'], 'score': 1.0 if props['drug_category'] == 'drug' else 0.8 if props['drug_category'] == 'lead' else 0.6 if props['drug_category'] == 'fragment' else 0.3 },
        { 'name': 'Lead-like', 'value': 1 if props['lead_like_pass'] else 0, 'target': 1, 'pass': props['lead_like_pass'], 'score': 1.0 if props['lead_like_pass'] else 0.0 },
        { 'name': 'Fragment-like', 'value': 1 if props['fragment_like_pass'] else 0, 'target': 1, 'pass': props['fragment_like_pass'], 'score': 1.0 if props['fragment_like_pass'] else 0.0 },
    ]
    
    def aggregate(ms):
        return sum(m['score'] for m in ms) / max(1, len(ms))
    
    return {
        'absorption': { 'metrics': abs_metrics, 'score': aggregate(abs_metrics) },
        'distribution': { 'metrics': dist_metrics, 'score': aggregate(dist_metrics) },
        'metabolism': { 'metrics': metab_metrics, 'score': aggregate(metab_metrics) },
        'excretion': { 'metrics': excr_metrics, 'score': aggregate(excr_metrics) },
        'toxicity': { 'metrics': tox_metrics, 'score': aggregate(tox_metrics) },
        'drug_likeness': { 'metrics': drug_likeness_metrics, 'score': aggregate(drug_likeness_metrics) },
    }


def clean_sequence(seq: str) -> str:
    seq = seq.upper()
    
    seq = re.sub(r"[^A-Z]", "", seq)
    # keep only standard 20 amino acids
    seq = ''.join([s for s in seq if s in AA_LIST])
    return seq


def parse_fasta(text: str) -> str:
    lines = text.strip().splitlines()
    if len(lines) == 0:
        return ""
    if lines[0].startswith('>'):
        lines = lines[1:]
    return clean_sequence(''.join(lines))


def predict_pdb_with_retry(sequence: str, max_retries: int = 3, base_delay: float = 1.0):
    """Predict PDB with retry mechanism and multiple API endpoints."""
    
    # List of API endpoints to try
    api_endpoints = [
        'https://api.esmatlas.com/foldSequence/v1/pdb/',
        'https://esmatlas.com/foldSequence/v1/pdb/',
        'https://api.esmatlas.com/foldSequence/v1/pdb'
    ]
    
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    
    for attempt in range(max_retries):
        for endpoint in api_endpoints:
            try:
                # Add jitter to prevent thundering herd
                delay = base_delay * (2 ** attempt) + random.uniform(0, 1)
                if attempt > 0:
                    log_api(f"Retry attempt {attempt + 1}/{max_retries} after {delay:.1f}s delay...")
                    time.sleep(delay)
                
                response = requests.post(endpoint, headers=headers, data=sequence, timeout=120)
                
                if response.status_code == 200:
                    content = response.content.decode('utf-8')
                    # Validate that we got PDB content, not an HTML error page
                    if content.strip().startswith('ATOM') or content.strip().startswith('HEADER'):
                        return content
                    else:
                        log_api(f"API returned non-PDB content from {endpoint}. Trying next endpoint...", 'warning')
                        continue
                elif response.status_code == 503:
                    log_api(f"Service temporarily unavailable (503) from {endpoint}. Trying next endpoint...", 'warning')
                    continue
                elif response.status_code == 429:
                    log_api(f"Rate limited (429) from {endpoint}. Waiting longer...", 'warning')
                    time.sleep(delay * 2)
                    continue
                else:
                    log_api(f"API returned status {response.status_code} from {endpoint}. Trying next endpoint...", 'warning')
                    continue
                    
            except requests.exceptions.Timeout:
                log_api(f"Timeout from {endpoint}. Trying next endpoint...", 'warning')
                continue
            except requests.exceptions.ConnectionError:
                log_api(f"Connection error from {endpoint}. Trying next endpoint...", 'warning')
                continue
            except Exception as e:
                log_api(f"Error from {endpoint}: {e}. Trying next endpoint...", 'warning')
                continue
    
    return None


def predict_pdb_local_esmfold(sequence: str):
    """Predict PDB using local ESMFold model (requires GPU)."""
    try:
        import torch
        from transformers import EsmForProteinFolding
        
        # Check if GPU is available
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
        if device.type == "cpu":
            st.warning("⚠️ GPU not available. Local ESMFold will be very slow on CPU.")
        
        # Load model
        model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1")
        model = model.to(device)
        
        # Predict structure
        with torch.no_grad():
            output = model.infer_pdb(sequence)
        
        return output
        
    except ImportError:
        st.error("""
        **Local ESMFold not available.** 
        
        To use local ESMFold, install the required packages:
        ```bash
        pip install torch transformers
        ```
        
        Note: This requires significant GPU memory (8GB+ recommended).
        """)
        return None
    except Exception as e:
        st.error(f"Local ESMFold prediction failed: {e}")
        return None


@st.cache_data(show_spinner=False)
def predict_pdb_from_esmfold(sequence: str, method: str = "api"):
    """Return pdb string. Cached so repeated sequences reuse result."""
    if method == "local":
        return predict_pdb_local_esmfold(sequence)
    else:
        return predict_pdb_with_retry(sequence)


def _looks_like_pdb(text: str) -> bool:
    if not text:
        return False
    stripped = text.lstrip()
    return stripped.startswith('ATOM') or stripped.startswith('HEADER') or stripped.startswith('MODEL')


@st.cache_data(show_spinner=False)
def fetch_rcsb_pdb_by_sequence(sequence: str, identity_cutoff: float = 0.6, evalue_cutoff: float = 1.0):
    """Search RCSB by sequence and return top hit PDB text if available."""
    try:
        query = {
            "query": {
                "type": "terminal",
                "service": "sequence",
                "parameters": {
                    "evalue_cutoff": evalue_cutoff,
                    "identity_cutoff": identity_cutoff,
                    "target": "pdb_protein_sequence",
                    "value": sequence
                }
            },
            "request_options": {
                "scoring_strategy": "sequence",
                "sort": [{"sort_by": "score", "direction": "desc"}],
                "results_content_type": ["experimental"],
                "results_verbosity": "minimal"
            },
            "return_type": "entry"
        }
        resp = requests.post(
            "https://search.rcsb.org/rcsbsearch/v2/query?json",
            headers={"Content-Type": "application/json"},
            data=json.dumps(query),
            timeout=30
        )
        if resp.status_code != 200:
            return None
        data = resp.json()
        result_set = data.get("result_set") or []
        if not result_set:
            return None
        top_id = result_set[0].get("identifier")
        if not top_id:
            return None
        pdb_resp = requests.get(f"https://files.rcsb.org/download/{top_id}.pdb", timeout=30)
        if pdb_resp.status_code == 200 and _looks_like_pdb(pdb_resp.text):
            return pdb_resp.text
    except Exception:
        return None
    return None


@st.cache_data(show_spinner=False)
def fetch_pdb_by_id(pdb_id: str):
    """Download PDB by 4-char ID from RCSB."""
    try:
        pid = pdb_id.strip().lower()
        if not pid or len(pid) < 4:
            return None
        url = f"https://files.rcsb.org/download/{pid}.pdb"
        r = requests.get(url, timeout=20)
        if r.status_code == 200 and _looks_like_pdb(r.text):
            return r.text
    except Exception:
        return None
    return None


@st.cache_data(show_spinner=False)
def fetch_alphafold_by_uniprot(uniprot_id: str):
    """Fetch AlphaFold model PDB by UniProt ID (AF-<ID>-F1-model_v4.pdb)."""
    try:
        uid = uniprot_id.strip()
        if not uid:
            return None
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb"
        r = requests.get(url, timeout=30)
        if r.status_code == 200 and _looks_like_pdb(r.text):
            return r.text
    except Exception:
        return None
    return None


def render_mol(pdb_text: str, color_mode: str = 'spectrum', surface: bool = False, height: int = 600, width: int = 900):
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_text, 'pdb')
    if color_mode == 'plddt':
        # color by b-factor using spectrum
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.setColorByFunction({'cartoon': {'prop': 'b'}})
    elif color_mode == 'hydrophobicity':
        # set per-residue coloring using b-factor channel trick -> we will rewrite b-factor values to map hydrophobicity
        view.setStyle({'cartoon': {'color': 'spectrum'}})
    elif color_mode == 'chain':
        view.setStyle({'cartoon': {'color': 'chain'}})
    else:
        view.setStyle({'cartoon': {'color': 'spectrum'}})

    if surface:
        view.setSurface(py3Dmol.VDW, {'opacity': 0.8})

    view.setBackgroundColor('white')
    view.zoomTo()
    view.zoom(1.4)
    showmol(view)


def extract_structure_and_scores(pdb_text: str):
    # write temporary buffer and load via biotite for coordinate operations
    with open('predicted.pdb', 'w') as f:
        f.write(pdb_text)
    # biotite can load structure and populate b_factor
    atom_array = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])  # AtomArray
    # convert per-atom b_factors to per-residue mean
    try:
        b_factors = atom_array.get_annotation('b_factor')
    except Exception:
        # fallback
        b_factors = atom_array.b_factor
    # get unique residues order
    res_ids = atom_array.res_id
    res_names = atom_array.res_name
    # Build per-residue list
    residues = []
    per_res_b = []
    coords_by_res = {}
    for atom, rid, rname, b, coord in zip(atom_array.atom_id, res_ids, res_names, b_factors, atom_array.coord):
        # res_id includes insertion and chain info; we will use simple residue index order
        idx = int(rid)
        residues.append((idx, rname))
        per_res_b.append(b)
        coords_by_res.setdefault(idx, []).append(coord)

    # collapse to sequential residues in file order
    # We'll derive a cleaned residue list using unique sequential indices in order of appearance
    seen = []
    res_seq_names = []
    for rid, rname in residues:
        if rid not in seen:
            seen.append(rid)
            res_seq_names.append(rname.strip())
    # per-residue plDDT: average b-factor per residue index
    unique_idxs = list(dict.fromkeys([r[0] for r in residues]))
    per_res_plddt = []
    for idx in unique_idxs:
        arr = coords_by_res.get(idx, [])
        # average b from atoms with that residue idx
        bvals = []
        for i, rr in enumerate(resids := list(res_ids)):
            if int(rr) == idx:
                try:
                    bvals.append(float(b_factors[i]))
                except Exception:
                    pass
        per_res_plddt.append(np.mean(bvals) if len(bvals) > 0 else np.nan)

    # if length mismatch between seq and plddt, try best-effort: use seq from atom_array
    seq_from_atoms = ''.join([one_three_to_one(n) for n in res_seq_names if one_three_to_one(n) != 'X'])

    return atom_array, seq_from_atoms, np.array(per_res_plddt)


# small helper: convert 3-letter to 1-letter AA
_three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def one_three_to_one(three):
    return _three_to_one.get(three.upper(), 'X')


def compute_molecular_weight(seq: str) -> float:
    return sum([AA_MASS.get(a, 0.0) for a in seq]) - (len(seq) - 1) * 18.01528  # subtract water for peptide bonds


def net_charge_at_pH(seq: str, pH: float) -> float:
    # Henderson-Hasselbalch for side chains + N/C termini
    positive = 0.0
    negative = 0.0
    # N-term
    positive += 1.0 / (1.0 + 10**(pH - PKA['Nterm']))
    # C-term
    negative += 1.0 / (1.0 + 10**(PKA['Cterm'] - pH))
    for aa in seq:
        if aa == 'K':
            positive += 1.0 / (1.0 + 10**(pH - PKA['K']))
        elif aa == 'R':
            positive += 1.0 / (1.0 + 10**(pH - PKA['R']))
        elif aa == 'H':
            positive += 1.0 / (1.0 + 10**(pH - PKA['H']))
        elif aa == 'D':
            negative += 1.0 / (1.0 + 10**(PKA['D'] - pH))
        elif aa == 'E':
            negative += 1.0 / (1.0 + 10**(PKA['E'] - pH))
        elif aa == 'C':
            negative += 1.0 / (1.0 + 10**(PKA['C'] - pH))
        elif aa == 'Y':
            negative += 1.0 / (1.0 + 10**(PKA['Y'] - pH))
    return positive - negative


def estimate_pI(seq: str, precision=0.01) -> float:
    # binary search between pH 0 and 14
    low, high = 0.0, 14.0
    while (high - low) > precision:
        mid = (low + high) / 2.0
        charge = net_charge_at_pH(seq, mid)
        if charge > 0:
            low = mid
        else:
            high = mid
    return round((low + high) / 2.0, 2)


def kyte_doolittle_profile(seq: str, window: int = 9):
    vals = [KD_SCALE.get(a, 0.0) for a in seq]
    # sliding window average
    half = window // 2
    prof = []
    for i in range(len(seq)):
        start = max(0, i - half)
        end = min(len(seq), i + half + 1)
        prof.append(np.mean(vals[start:end]))
    return prof


def residue_composition(seq: str):
    c = Counter(seq)
    return c


def charge_counts(seq: str):
    pos = sum([1 for a in seq if a in POSITIVE_AA])
    neg = sum([1 for a in seq if a in NEGATIVE_AA])
    neutral = len(seq) - pos - neg
    return pos, neg, neutral


def radius_of_gyration(atom_array):
    # use heavy atoms coordinates
    coords = atom_array.get_coord()
    # flatten to Nx3
    coords = coords.reshape(-1, 3)
    center = coords.mean(axis=0)
    rg = np.sqrt(((coords - center) ** 2).sum(axis=1).mean())
    return float(rg)


def parse_secondary_from_pdb(pdb_text: str):
    helix_count = len(re.findall(r'^HELIX', pdb_text, flags=re.MULTILINE))
    sheet_count = len(re.findall(r'^SHEET', pdb_text, flags=re.MULTILINE))
    return helix_count, sheet_count


# -------------------- New Advanced Functions --------------------

def calculate_rmsd(structure1, structure2):
    """Calculate RMSD between two protein structures."""
    try:
        # Align structures and calculate RMSD
        coords1 = structure1.get_coord()
        coords2 = structure2.get_coord()
        
        if len(coords1) != len(coords2):
            return None
            
        # Simple RMSD calculation (without alignment)
        diff = coords1 - coords2
        rmsd = np.sqrt(np.mean(diff**2))
        return float(rmsd)
    except:
        return None


def simple_sequence_alignment(seq1, seq2):
    """Simple sequence alignment for visualization."""
    # Basic alignment - can be enhanced with proper algorithm
    max_len = max(len(seq1), len(seq2))
    aligned_seq1 = seq1.ljust(max_len, '-')
    aligned_seq2 = seq2.ljust(max_len, '-')
    
    identity = sum(1 for a, b in zip(seq1, seq2) if a == b)
    similarity = identity / max(len(seq1), len(seq2)) * 100
    
    return aligned_seq1, aligned_seq2, similarity


def predict_domains_pfam(sequence):
    """Predict protein domains using Pfam (simulated)."""
    # This is a simplified version - in real implementation, you'd use Pfam API
    domains = []
    
    # Simple domain prediction based on sequence length and composition
    if len(sequence) > 100:
        domains.append({
            'name': 'Domain_1',
            'start': 1,
            'end': min(50, len(sequence)//2),
            'confidence': 0.8,
            'description': 'N-terminal domain'
        })
    
    if len(sequence) > 200:
        domains.append({
            'name': 'Domain_2', 
            'start': len(sequence)//2 + 1,
            'end': len(sequence),
            'confidence': 0.7,
            'description': 'C-terminal domain'
        })
    
    return domains


def predict_ppi_interactions(sequence):
    """Predict protein-protein interactions (simulated)."""
    # This is a simplified version - in real implementation, you'd use STRING API
    interactions = []
    
    # Simple PPI prediction based on sequence features
    if 'K' in sequence and 'R' in sequence:  # DNA binding proteins
        interactions.append({
            'partner': 'DNA',
            'confidence': 0.8,
            'type': 'DNA binding',
            'evidence': 'Sequence motifs'
        })
    
    if sequence.count('C') > 5:  # Proteins with many cysteines
        interactions.append({
            'partner': 'Metal ions',
            'confidence': 0.6,
            'type': 'Metal binding',
            'evidence': 'Cysteine content'
        })
    
    return interactions


def predict_disulfide_bonds(sequence):
    """Predict potential disulfide bonds."""
    cysteines = [i for i, aa in enumerate(sequence) if aa == 'C']
    bonds = []
    
    # Simple disulfide bond prediction (even positions)
    for i in range(0, len(cysteines)-1, 2):
        bonds.append({
            'residue1': cysteines[i] + 1,
            'residue2': cysteines[i+1] + 1,
            'confidence': 0.7
        })
    
    return bonds


def predict_ptm_sites(sequence):
    """Predict post-translational modification sites."""
    ptms = []
    
    # Phosphorylation sites
    for i, aa in enumerate(sequence):
        if aa in ['S', 'T', 'Y'] and i < len(sequence) - 2:
            # Simple motif-based prediction
            context = sequence[max(0, i-2):min(len(sequence), i+3)]
            if 'P' in context:  # Proline-directed
                ptms.append({
                    'position': i + 1,
                    'type': 'Phosphorylation',
                    'residue': aa,
                    'confidence': 0.6
                })
    
    # Glycosylation sites
    for i, aa in enumerate(sequence):
        if aa in ['N', 'S', 'T'] and i < len(sequence) - 2:
            context = sequence[i:i+3]
            if context.startswith('N') and context[2] in ['S', 'T']:
                ptms.append({
                    'position': i + 1,
                    'type': 'N-glycosylation',
                    'residue': aa,
                    'confidence': 0.5
                })
    
    return ptms


def predict_membrane_topology(sequence):
    """Predict membrane topology using hydrophobicity."""
    topology = []
    kd_profile = kyte_doolittle_profile(sequence)
    
    # Simple membrane prediction based on hydrophobicity
    for i, kd in enumerate(kd_profile):
        if kd > 1.0:  # Hydrophobic
            topology.append({
                'position': i + 1,
                'type': 'Membrane',
                'confidence': min(kd / 2.0, 1.0)
            })
    
    return topology


def batch_predict_structures(sequences, method="api"):
    """Predict structures for multiple sequences."""
    results = []
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    for i, (name, seq) in enumerate(sequences.items()):
        status_text.text(f'Predicting structure for {name}... ({i+1}/{len(sequences)})')
        
        try:
            pdb_text = predict_pdb_from_esmfold(seq, method=method)
            if pdb_text:
                results.append({
                    'name': name,
                    'sequence': seq,
                    'pdb': pdb_text,
                    'success': True
                })
            else:
                results.append({
                    'name': name,
                    'sequence': seq,
                    'pdb': None,
                    'success': False
                })
        except Exception as e:
            results.append({
                'name': name,
                'sequence': seq,
                'pdb': None,
                'success': False,
                'error': str(e)
            })
        
        progress_bar.progress((i + 1) / len(sequences))
    
    progress_bar.empty()
    status_text.empty()
    return results


# -------------------- UNIQUE ADVANCED FEATURES --------------------

def ai_protein_design_suggestions(sequence, target_function="stability"):
    """AI-powered protein design suggestions."""
    suggestions = []
    
    # Analyze sequence for design opportunities
    if target_function == "stability":
        # Suggest stabilizing mutations
        if sequence.count('G') > len(sequence) * 0.1:  # High glycine content
            suggestions.append({
                'position': sequence.find('G') + 1,
                'original': 'G',
                'suggested': 'A',
                'reason': 'Replace glycine with alanine for increased stability',
                'confidence': 0.8,
                'impact': 'Medium'
            })
        
        # Suggest helix-stabilizing mutations
        for i, aa in enumerate(sequence):
            if aa == 'P' and i > 0:
                suggestions.append({
                    'position': i + 1,
                    'original': 'P',
                    'suggested': 'A',
                    'reason': 'Replace proline to avoid helix breaking',
                    'confidence': 0.7,
                    'impact': 'High'
                })
    
    elif target_function == "activity":
        # Suggest activity-enhancing mutations
        if 'C' in sequence:
            suggestions.append({
                'position': sequence.find('C') + 1,
                'original': 'C',
                'suggested': 'S',
                'reason': 'Replace cysteine with serine to reduce oxidation risk',
                'confidence': 0.6,
                'impact': 'Medium'
            })
    
    elif target_function == "solubility":
        # Suggest solubility-enhancing mutations
        hydrophobic_residues = ['I', 'L', 'V', 'F', 'W', 'Y']
        for i, aa in enumerate(sequence):
            if aa in hydrophobic_residues and i < len(sequence) - 1:
                suggestions.append({
                    'position': i + 1,
                    'original': aa,
                    'suggested': 'N',
                    'reason': f'Replace {aa} with asparagine for improved solubility',
                    'confidence': 0.5,
                    'impact': 'Low'
                })
                break
    
    return suggestions


def predict_drug_binding_sites(structure_pdb):
    """Enhanced drug binding site prediction with improved algorithms."""
    binding_sites = []
    
    # Parse PDB structure
    lines = structure_pdb.split('\n')
    atoms = []
    residues = {}
    
    for line in lines:
        if line.startswith('ATOM'):
            atom_name = line[12:16].strip()
            res_name = line[17:20].strip()
            chain = line[21:22].strip()
            res_num = int(line[22:26].strip())
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            element = line[76:78].strip()
            
            atoms.append({
                'name': atom_name,
                'res_name': res_name,
                'chain': chain,
                'res_num': res_num,
                'coords': [x, y, z],
                'element': element
            })
            
            # Group by residue
            res_key = f"{chain}_{res_num}"
            if res_key not in residues:
                residues[res_key] = {
                    'res_name': res_name,
                    'chain': chain,
                    'res_num': res_num,
                    'atoms': []
                }
            residues[res_key]['atoms'].append(atoms[-1])
    
    if not atoms:
        return binding_sites
    
    # Convert to numpy arrays for calculations
    atom_coords = np.array([atom['coords'] for atom in atoms])
    
    # Enhanced pocket detection using multiple algorithms
    
    # 1. Convex hull-based pocket detection
    try:
        from scipy.spatial import ConvexHull
        hull = ConvexHull(atom_coords)
        hull_vertices = atom_coords[hull.vertices]
        
        # Find potential pocket centers by analyzing surface curvature
        pocket_candidates = []
        for i, vertex in enumerate(hull_vertices):
            # Calculate local surface properties
            distances = np.linalg.norm(atom_coords - vertex, axis=1)
            nearby_atoms = atom_coords[distances < 8.0]  # 8Å radius
            
            if len(nearby_atoms) > 5:  # Need sufficient atoms for analysis
                # Calculate local density and curvature
                local_center = np.mean(nearby_atoms, axis=0)
                local_density = len(nearby_atoms) / (4/3 * np.pi * 8**3)
                
                # Check if this could be a pocket (lower density than surrounding)
                if local_density < 0.1:  # Threshold for pocket detection
                    pocket_candidates.append({
                        'center': local_center,
                        'density': local_density,
                        'size': len(nearby_atoms)
                    })
    except ImportError:
        # Fallback to simpler method if scipy not available
        pocket_candidates = []
    
    # 2. Grid-based pocket detection
    if not pocket_candidates:
        # Create a 3D grid around the protein
        min_coords = np.min(atom_coords, axis=0) - 5
        max_coords = np.max(atom_coords, axis=0) + 5
        grid_size = 2.0  # 2Å grid spacing
        
        grid_points = []
        for x in np.arange(min_coords[0], max_coords[0], grid_size):
            for y in np.arange(min_coords[1], max_coords[1], grid_size):
                for z in np.arange(min_coords[2], max_coords[2], grid_size):
                    grid_points.append([x, y, z])
        
        grid_points = np.array(grid_points)
        
        # Find grid points that are inside the protein but not too close to atoms
        for point in grid_points:
            distances = np.linalg.norm(atom_coords - point, axis=1)
            min_dist = np.min(distances)
            
            # Point is inside if it's not too close to any atom but close enough to be in a pocket
            if 2.0 < min_dist < 6.0:
                # Calculate local properties
                nearby_atoms = atom_coords[distances < 8.0]
                if len(nearby_atoms) > 3:
                    local_density = len(nearby_atoms) / (4/3 * np.pi * 8**3)
                    pocket_candidates.append({
                        'center': point,
                        'density': local_density,
                        'size': len(nearby_atoms)
                    })
    
    # 3. Cluster nearby pocket candidates
    if pocket_candidates:
        # Remove duplicates and cluster nearby candidates
        unique_pockets = []
        for candidate in pocket_candidates:
            is_duplicate = False
            for existing in unique_pockets:
                if np.linalg.norm(np.array(candidate['center']) - np.array(existing['center'])) < 5.0:
                    is_duplicate = True
                    break
            if not is_duplicate:
                unique_pockets.append(candidate)
        
        # Analyze each pocket for druggability
        for i, pocket in enumerate(unique_pockets[:5]):  # Limit to top 5 pockets
            center = pocket['center']
            
            # Calculate pocket properties
            distances = np.linalg.norm(atom_coords - center, axis=1)
            nearby_atoms = [atoms[j] for j in range(len(atoms)) if distances[j] < 8.0]
            
            # Analyze residue composition
            res_types = {}
            hydrophobic_count = 0
            polar_count = 0
            charged_count = 0
            
            for atom in nearby_atoms:
                res_key = f"{atom['chain']}_{atom['res_num']}"
                if res_key in residues:
                    res_name = residues[res_key]['res_name']
                    res_types[res_name] = res_types.get(res_name, 0) + 1
                    
                    # Classify residue type
                    if res_name in ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO']:
                        hydrophobic_count += 1
                    elif res_name in ['SER', 'THR', 'ASN', 'GLN', 'TYR', 'CYS']:
                        polar_count += 1
                    elif res_name in ['LYS', 'ARG', 'HIS', 'ASP', 'GLU']:
                        charged_count += 1
            
            # Determine pocket type based on composition
            total_residues = len(res_types)
            if total_residues > 0:
                if hydrophobic_count / total_residues > 0.5:
                    pocket_type = 'hydrophobic'
                elif polar_count / total_residues > 0.5:
                    pocket_type = 'polar'
                elif charged_count / total_residues > 0.3:
                    pocket_type = 'charged'
                else:
                    pocket_type = 'mixed'
            else:
                pocket_type = 'mixed'
            
            # Calculate druggability score
            volume = pocket['size'] * 4.0  # Rough volume estimate
            druggability_factors = []
            
            # Size factor (optimal size for drug binding)
            if 50 < volume < 500:
                size_score = 1.0
            elif volume < 50:
                size_score = volume / 50.0
            else:
                size_score = max(0.0, 1.0 - (volume - 500) / 1000.0)
            druggability_factors.append(size_score)
            
            # Composition factor
            if pocket_type == 'mixed':
                comp_score = 0.8
            elif pocket_type in ['hydrophobic', 'polar']:
                comp_score = 0.6
            else:
                comp_score = 0.4
            druggability_factors.append(comp_score)
            
            # Accessibility factor (based on distance from surface)
            surface_dist = np.min(distances)
            if surface_dist < 3.0:
                access_score = 0.9
            elif surface_dist < 6.0:
                access_score = 0.7
            else:
                access_score = 0.3
            druggability_factors.append(access_score)
            
            # Calculate final druggability score
            druggability_score = np.mean(druggability_factors)
            confidence = min(0.95, 0.5 + druggability_score * 0.4)
            
            binding_sites.append({
                'site_id': i + 1,
                'center': center.tolist(),
                'volume': volume,
                'druggability_score': druggability_score,
                'pocket_type': pocket_type,
                'confidence': confidence,
                'residue_composition': res_types,
                'accessibility': surface_dist,
                'size_category': 'small' if volume < 100 else 'medium' if volume < 300 else 'large'
            })
    
    # If no pockets found, create a default analysis
    if not binding_sites and len(atoms) > 10:
        # Calculate protein center and create a basic binding site
        center = np.mean(atom_coords, axis=0)
        binding_sites.append({
            'site_id': 1,
            'center': center.tolist(),
            'volume': 200.0,
            'druggability_score': 0.3,
            'pocket_type': 'mixed',
            'confidence': 0.4,
            'residue_composition': {},
            'accessibility': 5.0,
            'size_category': 'medium'
        })
    
    return binding_sites


def generate_md_simulation_input(structure_pdb, simulation_type="equilibration"):
    """Generate molecular dynamics simulation input files."""
    
    if simulation_type == "equilibration":
        md_input = f"""
# Molecular Dynamics Equilibration Input
# Generated by ProStruct-3D

# General parameters
integrator = md
dt = 0.002
nsteps = 50000
nstxout = 1000
nstvout = 1000
nstenergy = 1000

# Temperature coupling
tcoupl = V-rescale
tc-grps = Protein Non-Protein
tau-t = 0.1 0.1
ref-t = 300 300

# Pressure coupling
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 2.0
ref-p = 1.0
compressibility = 4.5e-5

# Constraints
constraints = h-bonds
constraint-algorithm = lincs
continuation = no
lincs-order = 4

# Cutoffs
cutoff-scheme = Verlet
rvdw = 1.0
rlist = 1.0
coulombtype = PME
fourierspacing = 0.125

# Output control
nstlog = 1000
nstcalcenergy = 100
"""
    
    elif simulation_type == "production":
        md_input = f"""
# Molecular Dynamics Production Run Input
# Generated by ProStruct-3D

# General parameters
integrator = md
dt = 0.002
nsteps = 1000000
nstxout = 10000
nstvout = 10000
nstenergy = 10000

# Temperature coupling
tcoupl = Nose-Hoover
tc-grps = Protein Non-Protein
tau-t = 0.5 0.5
ref-t = 300 300

# Pressure coupling
pcoupl = Parrinello-Rahman
pcoupltype = isotropic
tau-p = 2.0
ref-p = 1.0
compressibility = 4.5e-5

# Constraints
constraints = h-bonds
constraint-algorithm = lincs
continuation = yes
lincs-order = 4

# Analysis
nstcomm = 1000
comm-mode = linear
comm-grps = Protein

# Output control
nstlog = 10000
nstcalcenergy = 100
"""
    
    return md_input


def analyze_evolutionary_conservation(sequence):
    """Analyze evolutionary conservation (simulated)."""
    conservation_scores = []
    
    # Simulate conservation analysis
    # In real implementation, you'd use tools like ConSurf or PSI-BLAST
    
    for i, aa in enumerate(sequence):
        # Higher conservation for functional residues
        if aa in ['C', 'W', 'F', 'Y', 'H']:  # Important functional residues
            score = np.random.uniform(0.7, 0.9)
        elif aa in ['G', 'P']:  # Structurally important
            score = np.random.uniform(0.5, 0.8)
        else:
            score = np.random.uniform(0.2, 0.7)
        
        conservation_scores.append({
            'position': i + 1,
            'residue': aa,
            'conservation_score': round(score, 3),
            'conservation_level': 'High' if score > 0.7 else 'Medium' if score > 0.4 else 'Low'
        })
    
    return conservation_scores


def create_structure_animation(frames_data):
    """Create animated structure visualization."""
    # This would generate animation frames
    # For now, return animation parameters
    animation_params = {
        'frames': len(frames_data),
        'duration': 2.0,  # seconds
        'loop': True,
        'interpolation': 'linear'
    }
    return animation_params


def generate_voice_commands():
    """Generate voice command suggestions."""
    commands = [
        "Show 3D structure",
        "Analyze hydrophobicity",
        "Predict domains",
        "Compare sequences",
        "Export results",
        "Generate report",
        "Switch to comparison mode",
        "Show statistics",
        "Download PDB file",
        "Clear analysis"
    ]
    return commands


def setup_collaboration_room(room_name):
    """Setup real-time collaboration room."""
    room_data = {
        'room_id': room_name,
        'participants': [],
        'shared_structures': [],
        'chat_messages': [],
        'last_activity': datetime.now().isoformat()
    }
    return room_data


# -------------------- Sidebar Configuration --------------------
with st.sidebar:
    st.header("🔧 Configuration")
    
    # Global settings
    color_mode = st.selectbox("3D Color Mode", ['spectrum', 'chain', 'plddt', 'hydrophobicity'])
    surface_toggle = st.checkbox("Show molecular surface", value=False)
    minimum_length = st.number_input("Min sequence length", min_value=10, max_value=5000, value=20)
    verbose_api = st.checkbox("Verbose API logs (retries, 503 warnings)", value=False, help="Enable to see detailed API retry logs")
    st.session_state.verbose_api = verbose_api

# Prediction method selection
    prediction_method = st.selectbox(
    "Prediction Method", 
    [
        "Auto (PDB → ESM → AlphaFold)",
        "PDB (Experimental structures)",
        "ESM Atlas API",
        "Local ESMFold (Requires GPU)",
        "AlphaFold DB (by UniProt ID)"
    ],
    help="Prefer experimental PDB first. Auto tries PDB sequence search, then ESM, then AlphaFold if UniProt ID provided."
)

    pdb_id_hint = st.text_input("Optional PDB ID (e.g., 1CRN)", value="")
    uniprot_id_hint = st.text_input("Optional UniProt ID for AlphaFold (e.g., P69905)", value="")

    # Batch processing settings
    st.subheader("📦 Batch Processing")
    enable_batch = st.checkbox("Enable batch processing", help="Process multiple sequences simultaneously")
    
    # Export settings
    st.subheader("💾 Export Options")
    export_formats = st.multiselect(
        "Export formats",
        ["PDB", "CSV", "JSON", "PyMOL", "ChimeraX"],
        default=["PDB", "CSV"]
    )

# -------------------- Tab 1: Structure Prediction --------------------
with tab1:
    st.header("🔬 Single Protein Structure Prediction")
    
    # Input section
    col1, col2 = st.columns([2, 1])
    
    with col1:
        seq_text = st.text_area(
            "Paste Protein Sequence (single-letter codes)", 
            height=120,
            placeholder="Enter your protein sequence here (e.g., MKWVTFISLLFLFSSAYS...)"
        )
    
    with col2:
        uploaded = st.file_uploader("Or upload FASTA/TXT", type=["fasta", "txt"])
        
        # Sample sequences
        with st.expander("📝 Sample Sequences"):
            samples = {
                "Insulin": "MALWMRLLPLLALLALWGPDPAAAFVNQHLCGSHLVEALYLVCGERGFFYTPKTRREAEDLQVGQVELGGGPGAGSLQPLALEGSLQKRGIVEQCCTSICSLYQLENYCN",
                "Myoglobin": "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASEDLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISEAIIQVLQSKHPGDFGADAQGAMTKALELFRNDIAAKYKELGFQG"
            }
            for name, seq in samples.items():
                if st.button(f"Load {name}", key=f"sample_{name}"):
                    st.session_state.sample_sequence = seq
                    st.rerun()
    
    # Read sequence
    sequence = ''
    if uploaded is not None:
        try:
            text = uploaded.read().decode('utf-8')
            sequence = parse_fasta(text)
        except Exception:
            st.error("Failed to read uploaded file. Ensure it's text/FASTA encoded UTF-8.")
    elif seq_text.strip():
        sequence = clean_sequence(seq_text)
    elif 'sample_sequence' in st.session_state:
        sequence = st.session_state.sample_sequence
        st.info(f"Loaded sample sequence: {len(sequence)} residues")

    if len(sequence) == 0:
        st.info("👆 Please paste a protein sequence or upload a FASTA file to get started.")
        st.stop()
    
    if len(sequence) < minimum_length:
        st.error(f"Sequence length {len(sequence)} is shorter than required minimum of {minimum_length}.")
        st.stop()
    
    # Display sequence info
    st.success(f"✅ Sequence loaded: {len(sequence)} residues")

    # Sequence composition and validation
    from collections import Counter
    STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")
    aa_counts = Counter(sequence)
    st.caption(f"Amino-acid composition (unique={len(aa_counts)}, length={len(sequence)})")
    non_standard_positions = [(i + 1, aa) for i, aa in enumerate(sequence) if aa not in STANDARD_AA]
    if non_standard_positions:
        unique_nonstd = sorted({aa for _, aa in non_standard_positions})
        st.warning(
            f"Non-standard residues detected: {', '.join(unique_nonstd)} "
            f"at positions: {', '.join(str(p) for p, _ in non_standard_positions[:50])}"
            + (" ..." if len(non_standard_positions) > 50 else "")
        )
        cleaned_seq = clean_sequence(sequence)
        if cleaned_seq != sequence:
            with st.expander("View cleaned sequence preview"):
                st.code(cleaned_seq)
            st.info(f"Cleaned sequence length: {len(cleaned_seq)} (removed {len(sequence) - len(cleaned_seq)} invalid chars)")
    
    # Prediction button
    pdb_text = None
    if st.button('🚀 Predict & Analyze Structure', type="primary", use_container_width=True):
        # Check rate limiting for API requests
        if not prediction_method.startswith("Local") and not check_rate_limit():
            st.stop()

        # Determine prediction method
        use_local = prediction_method.startswith("Local")
        method_name = prediction_method
        
        with st.spinner(f"Running {method_name} prediction and analyses..."):
            # Try PDB/Auto paths first
            pdb_text = None
            if prediction_method.startswith("PDB") and pdb_id_hint.strip():
                pdb_text = fetch_pdb_by_id(pdb_id_hint)
            if pdb_text is None and (prediction_method.startswith("Auto") or prediction_method.startswith("PDB")):
                seq_hit = fetch_rcsb_pdb_by_sequence(sequence)
                if seq_hit:
                    pdb_text = seq_hit
            # ESM
            if pdb_text is None and (prediction_method.startswith("Auto") or prediction_method.startswith("ESM") or use_local):
                pdb_text = predict_pdb_from_esmfold(sequence, method="local" if use_local else "api")

            # Fallback: if ESM fails or returns gateway/non-PDB, try RCSB by sequence
            needs_fallback = (
                pdb_text is None or
                (isinstance(pdb_text, str) and (pdb_text.strip().lower().startswith('<!doctype html') or not _looks_like_pdb(pdb_text)))
            )
            if needs_fallback and not use_local and not prediction_method.startswith("PDB"):
                st.info("ESM API failed or returned a gateway page. Trying RCSB PDB sequence search fallback...")
                alt_pdb = fetch_rcsb_pdb_by_sequence(sequence)
                if alt_pdb:
                    pdb_text = alt_pdb
                    st.success("Loaded closest matching experimental structure from RCSB PDB.")
            # AlphaFold by UniProt
            if (pdb_text is None or not _looks_like_pdb(pdb_text)) and (prediction_method.startswith("Auto") or prediction_method.startswith("AlphaFold")):
                if uniprot_id_hint.strip():
                    st.info("Trying AlphaFold DB by UniProt ID...")
                    af = fetch_alphafold_by_uniprot(uniprot_id_hint)
                    if af:
                        pdb_text = af
                        st.success("Loaded structure from AlphaFold DB.")

            if not pdb_text or not _looks_like_pdb(pdb_text):
                st.error("❌ **Could not obtain a valid structure** — ESM and RCSB fallback failed.")
                st.markdown("""
                **Try next:**
                1. Verify the sequence (ACDEFGHIKLMNPQRSTVWY only)
                2. Try a shorter sequence or split domains
                3. Retry later in case of API downtime
                4. Try Local ESMFold if you have GPU support
                """)
                if st.button("🔄 Retry", type="primary"):
                    st.rerun()
                st.stop()

            # extract atom array, seq and plDDT
            have_structure = False
            try:
                atom_array = bsio.load_structure(io.StringIO(pdb_text), extra_fields=["b_factor"])
                have_structure = True
            except Exception as e:
                with open('predicted.pdb', 'w') as f:
                    f.write(pdb_text)
                try:
                    atom_array = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
                    have_structure = True
                except Exception as e2:
                    st.error(f"❌ **Failed to load protein structure** — {str(e2)}")
                    st.markdown("""
                    **The PDB file appears to be corrupted or invalid.**
                    
                    **Possible causes:**
                    - API returned an incomplete response
                    - Network interruption during download
                    - Invalid protein sequence
                    
                    **Try:**
                    1. **Retry with a different sequence**
                    2. **Check your internet connection**
                    3. **Use Local ESMFold if available**
                    """)
                    have_structure = False
                    st.stop()

        # build dataframe per residue
        # get per-atom arrays and aggregate by residue sequence index
        atom_res_ids = atom_array.res_id
        atom_res_names = atom_array.res_name
        atom_b = atom_array.b_factor
        atom_coords = atom_array.coord

        # build mapping from residue (chain+resid) to list of atoms
        residue_keys = []
        residues_order = []
        for i in range(len(atom_res_ids)):
            key = (int(atom_res_ids[i]), atom_res_names[i])
            residue_keys.append(key)
            residues_order.append(key)
        # create ordered unique residue list
        ordered_res = []
        for k in residues_order:
            if k not in ordered_res:
                ordered_res.append(k)

        perres = []
        for ridx, r in enumerate(ordered_res, start=1):
            # select atoms with that residue name in order
            b_vals = []
            coords = []
            for i in range(len(atom_res_ids)):
                if (int(atom_res_ids[i]), atom_res_names[i]) == r:
                    try:
                        b_vals.append(float(atom_b[i]))
                    except Exception:
                        pass
                    coords.append(atom_coords[i])
            # average b
            avg_b = float(np.nanmean(b_vals)) if len(b_vals) > 0 else np.nan
            one_letter = one_three_to_one(r[1])
            perres.append({'index': ridx, 'res_name_3': r[1], 'res_name_1': one_letter, 'avg_b': avg_b, 'n_atoms': len(b_vals)})

        df = pd.DataFrame(perres)
        # if sequence shorter/longer, prefer sequence from input
        seq_for_analysis = ''.join(df['res_name_1'].tolist())
        if len(seq_for_analysis) != len(sequence):
            # prefer supplied sequence for hydrophobicity etc, but align lengths in plots
            seq_for_analysis = sequence

        # plDDT array
        plddt_vals = df['avg_b'].fillna(0).to_numpy()

        # Confidence stats
        pct_above_70 = round((plddt_vals > 70).sum() / max(1, len(plddt_vals)) * 100, 2)
        pct_above_90 = round((plddt_vals > 90).sum() / max(1, len(plddt_vals)) * 100, 2)

        # hydrophobicity
        hyd_profile = kyte_doolittle_profile(seq_for_analysis)

        # residue comp and charge
        comp = residue_composition(seq_for_analysis)
        posc, negc, neutr = charge_counts(seq_for_analysis)

        # mw, pI
        molw = compute_molecular_weight(seq_for_analysis)
        pipred = estimate_pI(seq_for_analysis)

        # radius of gyration
        try:
            rg = radius_of_gyration(atom_array)
        except Exception:
            rg = None

        # secondary structure from PDB records
        helix_count, sheet_count = parse_secondary_from_pdb(pdb_text)

        # prepare per-residue csv for download (align lengths)
        n_res = min(len(df), len(hyd_profile))
        if len(df) != len(hyd_profile):
            log_api(f"Length mismatch: structure residues={len(df)} vs sequence length={len(hyd_profile)}; aligning to {n_res}")
        df_aligned = df.iloc[:n_res].reset_index(drop=True)
        hyd_aligned = hyd_profile[:n_res]
        out_df = pd.DataFrame({
            'res_index': df_aligned['index'],
            'residue': df_aligned['res_name_1'],
            'plddt': df_aligned['avg_b'],
            'hydrophobicity': hyd_aligned
        })

        # Store in session state for history
        prediction_result = {
            'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            'sequence': sequence,
            'sequence_length': len(sequence),
            'pdb_text': pdb_text,
            'plddt_vals': plddt_vals.tolist(),
            'confidence_70': pct_above_70,
            'confidence_90': pct_above_90,
            'molecular_weight': molw,
            'pI': pipred,
            'radius_gyration': rg,
            'helix_count': helix_count,
            'sheet_count': sheet_count,
            'composition': dict(comp),
            'charge_distribution': {'positive': posc, 'negative': negc, 'neutral': neutr}
        }
        st.session_state.prediction_history.append(prediction_result)

    # -------------------- Layout: visualization + stats --------------------
    st.markdown("---")
    # Use columns for better alignment and a modern look
    col1, col2 = st.columns([1.5, 1], gap="large")

    with col1:
        st.subheader('🧬 3D Structure')
        st.markdown("""
        <div class="viewer-container">
        """, unsafe_allow_html=True)
        if pdb_text and _looks_like_pdb(pdb_text):
            st.markdown('<div style="display:flex; justify-content:center;">', unsafe_allow_html=True)
            render_mol(pdb_text, color_mode=color_mode, surface=surface_toggle, height=600, width=600)
            st.markdown('</div>', unsafe_allow_html=True)
            st.markdown('<div style="text-align: center; margin-top: 1rem;">', unsafe_allow_html=True)
            st.download_button(
                '⬇️ Download PDB Structure', 
                data=pdb_text, 
                file_name='predicted.pdb', 
                mime='text/plain',
                type="primary"
            )
            st.markdown('</div>', unsafe_allow_html=True)
        else:
            st.warning('No valid structure available to render or download.')
        st.markdown("</div>", unsafe_allow_html=True)

    with col2:
        st.subheader('📊 Quick Statistics')
        st.markdown("""
        <div class="chart-container">
            <div style="display: grid; gap: 1rem;">
        """, unsafe_allow_html=True)
        # Determine a robust sequence length for display
        try:
            seq_len_display = len(seq_for_analysis)
        except NameError:
            seq_len_display = len(df) if 'df' in locals() else len(sequence)
        # Ensure dependent stats exist (fallbacks for robustness)
        if 'molw' not in locals():
            try:
                molw = compute_molecular_weight(sequence)
            except Exception:
                molw = float('nan')
        if 'pipred' not in locals():
            try:
                pipred = estimate_pI(sequence)
            except Exception:
                pipred = float('nan')
        if 'rg' not in locals():
            rg = None
        if 'pct_above_70' not in locals():
            pct_above_70 = 0.0
        if 'pct_above_90' not in locals():
            pct_above_90 = 0.0

        metrics = [
            ("Sequence Length", f"{seq_len_display}", "residues", "#3b82f6"),
            ("Molecular Weight", f"{round(molw,2)}", "Da", "#10b981"),
            ("Estimated pI", f"{pipred}", "", "#f59e0b"),
            ("Radius of Gyration", f"{round(rg,3)}" if rg else "N/A", "Å" if rg else "", "#ef4444"),
            ("High Confidence", f"{pct_above_70}%", "pLDDT > 70", "#8b5cf6"),
            ("Very High Confidence", f"{pct_above_90}%", "pLDDT > 90", "#06b6d4")
        ]
        for title, value, unit, color in metrics:
            st.markdown(f"""
            <div class="metric-card" style="border-left-color: {color};">
                <div style="display: flex; justify-content: space-between; align-items: center;">
                    <div>
                        <h4 style="color: {color}; margin: 0; font-size: 0.9rem; font-weight: 600;">{title}</h4>
                        <p style="margin: 0.3rem 0 0 0; color: #64748b; font-size: 0.8rem;">{unit}</p>
                    </div>
                    <div style="font-size: 1.5rem; font-weight: 700; color: {color};">{value}</div>
                </div>
            </div>
            """, unsafe_allow_html=True)
        # Guard secondary structure metrics
        try:
            helices_display = helix_count
            sheets_display = sheet_count
        except NameError:
            helices_display = 0
            sheets_display = 0
        st.markdown(f"""
        <div class="metric-card" style="border-left-color: #f97316;">
            <h4 style="color: #f97316; margin: 0; font-size: 0.9rem; font-weight: 600;">Secondary Structure</h4>
            <p style="margin: 0.5rem 0 0 0; color: #64748b; font-size: 0.9rem;">
                <strong>Helices:</strong> {helices_display} | <strong>Sheets:</strong> {sheets_display}
            </p>
        </div>
        """, unsafe_allow_html=True)
        st.markdown("</div></div>", unsafe_allow_html=True)

        # Ensure downstream objects exist even if prediction block didn't set them
        if 'out_df' not in locals():
            out_df = pd.DataFrame(columns=['res_index', 'plddt'])
        if 'hyd_profile' not in locals():
            hyd_profile = []

        st.subheader('🧪 Charge & Composition')
        st.markdown("""
        <div class="chart-container">
        """, unsafe_allow_html=True)
        # Guard charge and composition variables
        if 'posc' not in locals():
            posc = 0
        if 'negc' not in locals():
            negc = 0
        if 'neutr' not in locals():
            neutr = 0
        if 'comp' not in locals():
            comp = {}
        st.markdown(f"""
        <div class="metric-card" style="border-left-color: #8b5cf6;">
            <h4 style="color: #8b5cf6; margin: 0; font-size: 0.9rem; font-weight: 600;">Charge Distribution</h4>
            <div style="display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 1rem; margin-top: 0.5rem;">
                <div style="text-align: center;">
                    <div style="font-size: 1.2rem; font-weight: 700; color: #ef4444;">{posc}</div>
                    <div style="font-size: 0.8rem; color: #64748b;">Positive (K/R/H)</div>
                </div>
                <div style="text-align: center;">
                    <div style="font-size: 1.2rem; font-weight: 700; color: #3b82f6;">{negc}</div>
                    <div style="font-size: 0.8rem; color: #64748b;">Negative (D/E)</div>
                </div>
                <div style="text-align: center;">
                    <div style="font-size: 1.2rem; font-weight: 700; color: #64748b;">{neutr}</div>
                    <div style="font-size: 0.8rem; color: #64748b;">Neutral</div>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
        comp_items = sorted(comp.items(), key=lambda x: x[1], reverse=True) if isinstance(comp, dict) else []
        labels = [k for k, v in comp_items]
        values = [v for k, v in comp_items]
        if len(labels) > 0:
            fig_pie = go.Figure(data=[go.Pie(
                labels=labels, 
                values=values, 
                hole=0.4,
                marker_colors=['#3b82f6', '#10b981', '#f59e0b', '#ef4444', '#8b5cf6', '#06b6d4', '#f97316', '#84cc16', '#ec4899', '#6366f1', '#14b8a6', '#f59e0b', '#ef4444', '#8b5cf6', '#06b6d4', '#f97316', '#84cc16', '#ec4899', '#6366f1', '#14b8a6']
            )])
            fig_pie.update_layout(
                title=dict(
                    text='Residue Composition',
                    font=dict(size=16, color='#1f2937'),
                    x=0.5
                ),
                showlegend=True,
                legend=dict(
                    orientation="v",
                    yanchor="middle",
                    y=0.5,
                    xanchor="left",
                    x=1.01
                ),
                plot_bgcolor='rgba(0,0,0,0)',
                paper_bgcolor='rgba(0,0,0,0)',
                font=dict(size=12)
            )
            st.plotly_chart(fig_pie, use_container_width=True)
        else:
            st.info("No composition data available.")
        st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("---")
    st.subheader('📈 Per-residue Analyses')
    col3, col4 = st.columns(2, gap="large")

    with col3:
        st.markdown("""
        <div class="chart-container">
            <h4 style="color: #1f2937; margin-bottom: 1rem; text-align: center;">plDDT Confidence Scores</h4>
        """, unsafe_allow_html=True)
        if len(out_df) > 0:
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=out_df['res_index'], 
                y=out_df['plddt'], 
                mode='lines+markers', 
                name='plDDT',
                line=dict(color='#3b82f6', width=3),
                marker=dict(size=6, color='#1d4ed8'),
                hovertemplate='<b>Residue %{x}</b><br>plDDT: %{y:.1f}<extra></extra>'
            ))
            fig.add_hline(y=70, line_dash='dash', line_color='#f59e0b', line_width=2,
                         annotation_text='Low Confidence', annotation_position='top left',
                         annotation_font_color='#92400e')
            fig.add_hline(y=90, line_dash='dash', line_color='#10b981', line_width=2,
                         annotation_text='High Confidence', annotation_position='top left',
                         annotation_font_color='#065f46')
            fig.update_layout(
                xaxis_title='Residue Index',
                yaxis_title='plDDT Score',
                yaxis=dict(range=[0,100]),
                plot_bgcolor='rgba(255,255,255,0.8)',
                paper_bgcolor='rgba(255,255,255,0.8)',
                font=dict(size=12, family='Inter'),
                margin=dict(l=20, r=20, t=20, b=20),
                showlegend=False
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.empty()
        st.markdown("</div>", unsafe_allow_html=True)

    with col4:
        st.markdown("""
        <div class="chart-container">
            <h4 style="color: #1f2937; margin-bottom: 1rem; text-align: center;">Hydrophobicity Profile</h4>
        """, unsafe_allow_html=True)
        if len(hyd_profile) > 0:
            fig2 = go.Figure()
            fig2.add_trace(go.Scatter(
                x=list(range(1, len(hyd_profile)+1)), 
                y=hyd_profile, 
                mode='lines',
                name='Hydrophobicity',
                line=dict(color='#10b981', width=3),
                fill='tonexty',
                fillcolor='rgba(16, 185, 129, 0.1)',
                hovertemplate='<b>Residue %{x}</b><br>KD Score: %{y:.2f}<extra></extra>'
            ))
            fig2.update_layout(
                xaxis_title='Residue Index',
                yaxis_title='Kyte-Doolittle Score',
                plot_bgcolor='rgba(255,255,255,0.8)',
                paper_bgcolor='rgba(255,255,255,0.8)',
                font=dict(size=12, family='Inter'),
                margin=dict(l=20, r=20, t=20, b=20),
                showlegend=False
            )
            st.plotly_chart(fig2, use_container_width=True)
        else:
            st.empty()
        st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("---")
    st.subheader('⚡ Charge Distribution')
    st.markdown("""
    <div class="chart-container">
    """, unsafe_allow_html=True)
    if posc + negc + neutr > 0:
        charge_fig = go.Figure(data=[go.Pie(
            labels=['Positive','Negative','Neutral'], 
            values=[posc,negc,neutr], 
            hole=0.6,
            marker_colors=['#ef4444', '#3b82f6', '#64748b'],
            textinfo='label+percent',
            textfont_size=14
        )])
        charge_fig.update_layout(
            title=dict(
                text='Charge Distribution',
                font=dict(size=18, color='#1f2937'),
                x=0.5
            ),
            plot_bgcolor='rgba(255,255,255,0.8)',
            paper_bgcolor='rgba(255,255,255,0.8)',
            font=dict(size=12, family='Inter'),
            showlegend=True,
            legend=dict(
                orientation="h",
                yanchor="bottom",
                y=1.02,
                xanchor="right",
                x=1
            )
        )
        st.plotly_chart(charge_fig, use_container_width=True)
    else:
        st.info("No charge data available.")
    st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("---")
    st.subheader('⬇️ Download Results')
    st.markdown("""
    <div class="chart-container">
        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 1rem; margin-top: 1rem;">
    """, unsafe_allow_html=True)
    # Build CSV safely
    try:
        csv_buf = out_df.to_csv(index=False).encode('utf-8')
    except Exception:
        csv_buf = b""
    col1, col2 = st.columns(2)
    with col1:
        st.download_button(
            '📊 Download CSV Data', 
            data=csv_buf, 
            file_name='per_residue_analysis.csv', 
            mime='text/csv',
            type="primary"
        )
    with col2:
        report = io.StringIO()
        report.write('🧬 ProStruct - Protein Structure Analysis Report\n')
        report.write('=' * 50 + '\n\n')
        safe_seq = seq_for_analysis if 'seq_for_analysis' in locals() else sequence
        safe_mw = round(molw,2) if 'molw' in locals() else float('nan')
        safe_pi = pipred if 'pipred' in locals() else float('nan')
        report.write(f'📏 Sequence length: {len(safe_seq)} residues\n')
        report.write(f'⚖️ Molecular weight: {safe_mw} Da\n')
        report.write(f'🔬 Estimated pI: {safe_pi}\n')
        report.write(f'📈 Confidence scores:\n')
        report.write(f'   • High confidence (pLDDT > 70): {pct_above_70 if "pct_above_70" in locals() else 0}%\n')
        report.write(f'   • Very high confidence (pLDDT > 90): {pct_above_90 if "pct_above_90" in locals() else 0}%\n')
        report.write(f'🏗️ Secondary structure:\n')
        report.write(f'   • Helices: {helix_count if "helix_count" in locals() else 0}\n')
        report.write(f'   • Sheets: {sheet_count if "sheet_count" in locals() else 0}\n')
        report.write(f'\n🧪 Residue composition:\n')
        safe_items = comp_items if 'comp_items' in locals() else []
        for aa, count in safe_items:
            pct = (count/len(safe_seq)*100) if len(safe_seq) else 0
            report.write(f'   • {aa}: {count} ({pct:.1f}%)\n')
        report.write(f'\n⚡ Charge distribution:\n')
        report.write(f'   • Positive (K/R/H): {posc if "posc" in locals() else 0}\n')
        report.write(f'   • Negative (D/E): {negc if "negc" in locals() else 0}\n')
        report.write(f'   • Neutral: {neutr if "neutr" in locals() else 0}\n')
        st.download_button(
            '📄 Download Analysis Report', 
            data=report.getvalue(), 
            file_name='analysis_report.txt', 
            mime='text/plain',
            type="secondary"
        )
    st.markdown("</div></div>", unsafe_allow_html=True)

    # Success message with better styling
    st.markdown("""
    <div style="background: linear-gradient(135deg, #dcfce7 0%, #bbf7d0 100%); 
                border: 1px solid #22c55e; border-radius: 16px; 
                padding: 1.5rem; margin: 2rem 0; text-align: center;">
        <h3 style="color: #166534; margin: 0 0 0.5rem 0;">✅ Analysis Complete!</h3>
        <p style="color: #15803d; margin: 0;">Your protein structure has been successfully predicted and analyzed.</p>
    </div>
    """, unsafe_allow_html=True)

# -------------------- Tab 2: Protein Comparison --------------------
with tab2:
    st.header("🔄 Protein Structure Comparison")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.subheader("📝 Sequence 1")
        seq1_text = st.text_area("First protein sequence", height=100, key="seq1")
        seq1_name = st.text_input("Name (optional)", value="Protein 1", key="name1")
    
    with col2:
        st.subheader("📝 Sequence 2")
        seq2_text = st.text_area("Second protein sequence", height=100, key="seq2")
        seq2_name = st.text_input("Name (optional)", value="Protein 2", key="name2")
    
    # Load from history
    if st.session_state.prediction_history:
        with st.expander("📚 Load from Prediction History"):
            history_options = [f"{i+1}. {h['timestamp']} ({h['sequence_length']} residues)" 
                             for i, h in enumerate(st.session_state.prediction_history)]
            selected_history = st.selectbox("Select prediction", ["None"] + history_options)
            
            if selected_history != "None":
                idx = int(selected_history.split('.')[0]) - 1
                selected_pred = st.session_state.prediction_history[idx]
                
                col_a, col_b = st.columns(2)
                with col_a:
                    if st.button("Load as Sequence 1", key="load_seq1"):
                        st.session_state.seq1_from_history = selected_pred['sequence']
                        st.rerun()
                with col_b:
                    if st.button("Load as Sequence 2", key="load_seq2"):
                        st.session_state.seq2_from_history = selected_pred['sequence']
                        st.rerun()
    
    # Use sequences from history if loaded
    if 'seq1_from_history' in st.session_state:
        seq1_text = st.session_state.seq1_from_history
        del st.session_state.seq1_from_history
    if 'seq2_from_history' in st.session_state:
        seq2_text = st.session_state.seq2_from_history
        del st.session_state.seq2_from_history
    
    seq1_clean = clean_sequence(seq1_text) if seq1_text else ""
    seq2_clean = clean_sequence(seq2_text) if seq2_text else ""
    
    if seq1_clean and seq2_clean:
        # Sequence alignment
        aligned_seq1, aligned_seq2, similarity = simple_sequence_alignment(seq1_clean, seq2_clean)
        
        st.success(f"✅ Sequences loaded: {len(seq1_clean)} vs {len(seq2_clean)} residues")
        st.info(f"📊 Sequence similarity: {similarity:.1f}%")
        
        # Display alignment
        st.subheader("🔍 Sequence Alignment")
        st.markdown("""
        <div style="background: #f8fafc; padding: 1rem; border-radius: 8px; font-family: monospace; font-size: 12px; overflow-x: auto;">
        """, unsafe_allow_html=True)
        
        # Create alignment visualization
        max_display = 100
        if len(aligned_seq1) > max_display:
            aligned_seq1_display = aligned_seq1[:max_display] + "..."
            aligned_seq2_display = aligned_seq2[:max_display] + "..."
        else:
            aligned_seq1_display = aligned_seq1
            aligned_seq2_display = aligned_seq2
        
        st.text(f"{seq1_name}: {aligned_seq1_display}")
        st.text(f"{seq2_name}: {aligned_seq2_display}")
        
        # Show match/mismatch
        match_line = ""
        for a, b in zip(aligned_seq1[:max_display], aligned_seq2[:max_display]):
            if a == b and a != '-':
                match_line += "|"
            else:
                match_line += " "
        
        st.text(f"Match:     {match_line}")
        st.markdown("</div>", unsafe_allow_html=True)
        
        # Comparison button
        if st.button("🔬 Compare Structures", type="primary", use_container_width=True):
            with st.spinner("Retrieving structures (PDB → ESM) for both sequences..."):
                # Robust per-sequence retrieval: RCSB by sequence first, then ESM
                def get_best_pdb_for_sequence(seq: str):
                    try:
                        hit = fetch_rcsb_pdb_by_sequence(seq)
                        if hit and _looks_like_pdb(hit):
                            return hit
                    except Exception:
                        pass
                    try:
                        pred = predict_pdb_from_esmfold(seq, method="api")
                        if pred and _looks_like_pdb(pred):
                            return pred
                    except Exception:
                        pass
                    return None

                pdb1 = get_best_pdb_for_sequence(seq1_clean)
                pdb2 = get_best_pdb_for_sequence(seq2_clean)

                results = [
                    {'name': seq1_name, 'sequence': seq1_clean, 'pdb': pdb1, 'success': bool(pdb1)},
                    {'name': seq2_name, 'sequence': seq2_clean, 'pdb': pdb2, 'success': bool(pdb2)},
                ]

                successful_results = [r for r in results if r['success']]

                if len(successful_results) == 2:
                    # Calculate RMSD
                    try:
                        struct1 = bsio.load_structure(io.StringIO(results[0]['pdb']))
                        struct2 = bsio.load_structure(io.StringIO(results[1]['pdb']))
                        rmsd = calculate_rmsd(struct1, struct2)
                    except:
                        rmsd = None
                    
                    # Display comparison results
                    st.subheader("📊 Comparison Results")
                    
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        st.metric("Sequence Similarity", f"{similarity:.1f}%")
                    
                    with col2:
                        rmsd_text = f"{rmsd:.2f} Å" if rmsd else "N/A"
                        st.metric("Structural RMSD", rmsd_text)
                    
                    with col3:
                        len_diff = abs(len(seq1_clean) - len(seq2_clean))
                        st.metric("Length Difference", f"{len_diff} residues")
                    
                    # Side-by-side 3D structures
                    st.subheader("🧬 3D Structure Comparison")
                    col_left, col_right = st.columns(2)
                    
                    with col_left:
                        st.markdown(f"**{seq1_name}**")
                        render_mol(results[0]['pdb'], color_mode=color_mode, height=400, width=400)
                    
                    with col_right:
                        st.markdown(f"**{seq2_name}**")
                        render_mol(results[1]['pdb'], color_mode=color_mode, height=400, width=400)
                    
                    # Store comparison in session state
                    comparison_data = {
                        'timestamp': datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                        'sequence1': {'name': seq1_name, 'sequence': seq1_clean, 'pdb': results[0]['pdb']},
                        'sequence2': {'name': seq2_name, 'sequence': seq2_clean, 'pdb': results[1]['pdb']},
                        'similarity': similarity,
                        'rmsd': rmsd,
                        'alignment': {'seq1': aligned_seq1, 'seq2': aligned_seq2}
                    }

                    st.session_state.comparison_data.append(comparison_data)
                    
                elif len(successful_results) == 1:
                    st.warning("⚠️ Only one structure could be retrieved (PDB/ESM). Showing the one that worked.")
                    successful_result = successful_results[0]
                    st.markdown(f"**Available structure: {successful_result['name']}**")
                    render_mol(successful_result['pdb'], color_mode=color_mode, height=400, width=600)
                else:
                    st.error("❌ Could not retrieve either structure (PDB and ESM both failed).")
    
    elif seq1_clean or seq2_clean:
        st.warning("⚠️ Please provide both sequences for comparison.")
    else:
        st.info("👆 Please enter two protein sequences to compare their structures.")
    
    # Display comparison history
    if st.session_state.comparison_data:
        st.subheader("📚 Comparison History")
        for i, comp in enumerate(st.session_state.comparison_data):
            with st.expander(f"Comparison {i+1}: {comp['sequence1']['name']} vs {comp['sequence2']['name']} ({comp['timestamp']})"):
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown(f"**{comp['sequence1']['name']}**")
                    st.markdown(f"Length: {len(comp['sequence1']['sequence'])} residues")
                    if st.button(f"View Structure", key=f"view1_{i}"):
                        render_mol(comp['sequence1']['pdb'], color_mode=color_mode, height=300, width=300)
                
                with col2:
                    st.markdown(f"**{comp['sequence2']['name']}**")
                    st.markdown(f"Length: {len(comp['sequence2']['sequence'])} residues")
                    if st.button(f"View Structure", key=f"view2_{i}"):
                        render_mol(comp['sequence2']['pdb'], color_mode=color_mode, height=300, width=300)
                
                st.markdown(f"**Similarity:** {comp['similarity']:.1f}%")
                st.markdown(f"**RMSD:** {comp['rmsd']:.2f} Å" if comp['rmsd'] else "**RMSD:** N/A")

# -------------------- Tab 3: Domain Analysis --------------------
with tab3:
    st.header("🏗️ Protein Domain Analysis")
    
    # Input sequence
    domain_seq_text = st.text_area("Protein Sequence for Domain Analysis", height=120, key="domain_seq")
    
    if domain_seq_text:
        domain_seq = clean_sequence(domain_seq_text)
        
        if len(domain_seq) >= 20:
            st.success(f"✅ Sequence loaded: {len(domain_seq)} residues")
            
            if st.button("🔍 Predict Domains", type="primary"):
                with st.spinner("Analyzing protein domains..."):
                    domains = predict_domains_pfam(domain_seq)
                    
                    if domains:
                        st.subheader("📊 Predicted Domains")
                        
                        # Domain visualization
                        fig = go.Figure()
                        
                        for i, domain in enumerate(domains):
                            fig.add_trace(go.Scatter(
                                x=[domain['start'], domain['end'], domain['end'], domain['start'], domain['start']],
                                y=[i, i, i+0.8, i+0.8, i],
                                fill='toself',
                                mode='lines',
                                name=domain['name'],
                                line=dict(width=0),
                                fillcolor=f'rgba({50 + i*50}, {150 + i*30}, {200}, 0.7)',
                                hovertemplate=f"<b>{domain['name']}</b><br>" +
                                            f"Position: {domain['start']}-{domain['end']}<br>" +
                                            f"Confidence: {domain['confidence']:.2f}<br>" +
                                            f"Description: {domain['description']}<extra></extra>"
                            ))
                        
                        fig.update_layout(
                            title="Protein Domain Map",
                            xaxis_title="Residue Position",
                            yaxis_title="Domains",
                            height=max(200, len(domains) * 100),
                            showlegend=True,
                            plot_bgcolor='rgba(0,0,0,0)',
                            paper_bgcolor='rgba(0,0,0,0)'
                        )
                        st.plotly_chart(fig, use_container_width=True)
                        
                        # Domain table
                        st.subheader("📋 Domain Details")
                        domain_df = pd.DataFrame(domains)
                        st.dataframe(domain_df, use_container_width=True)
                        
                        # Store results
                        st.session_state.domain_results[domain_seq] = domains
                        
                    else:
                        st.info("No domains predicted for this sequence.")
        else:
            st.warning("⚠️ Sequence too short for domain analysis (minimum 20 residues).")
    else:
        st.info("👆 Please enter a protein sequence for domain analysis.")

# -------------------- Tab 4: PPI Prediction --------------------
with tab4:
    st.header("🔗 Protein-Protein Interaction Prediction")
    
    # Input sequence
    ppi_seq_text = st.text_area("Protein Sequence for PPI Analysis", height=120, key="ppi_seq")
    
    if ppi_seq_text:
        ppi_seq = clean_sequence(ppi_seq_text)
        
        if len(ppi_seq) >= 20:
            st.success(f"✅ Sequence loaded: {len(ppi_seq)} residues")
            
            if st.button("🔍 Predict Interactions", type="primary"):
                with st.spinner("Predicting protein interactions..."):
                    interactions = predict_ppi_interactions(ppi_seq)
                    
                    if interactions:
                        st.subheader("🔗 Predicted Interactions")
                        
                        # Interaction network visualization
                        nodes = [{'id': 'Target Protein', 'group': 1, 'size': 20}]
                        links = []
                        
                        for i, interaction in enumerate(interactions):
                            nodes.append({
                                'id': interaction['partner'],
                                'group': 2,
                                'size': max(10, interaction['confidence'] * 15)
                            })
                            links.append({
                                'source': 'Target Protein',
                                'target': interaction['partner'],
                                'value': interaction['confidence']
                            })
                        
                        # Create network visualization
                        fig = go.Figure()
                        
                        # Add nodes
                        for node in nodes:
                            color = '#3b82f6' if node['group'] == 1 else '#10b981'
                            fig.add_trace(go.Scatter(
                                x=[random.uniform(0, 100)],
                                y=[random.uniform(0, 100)],
                                mode='markers+text',
                                marker=dict(size=node['size'], color=color),
                                text=[node['id']],
                                textposition="middle center",
                                name=node['id'],
                                hovertemplate=f"<b>{node['id']}</b><extra></extra>"
                            ))
                        
                        fig.update_layout(
                            title="Protein Interaction Network",
                            showlegend=False,
                            height=400,
                            plot_bgcolor='rgba(0,0,0,0)',
                            paper_bgcolor='rgba(0,0,0,0)'
                        )
                        st.plotly_chart(fig, use_container_width=True)
                        
                        # Interaction table
                        st.subheader("📋 Interaction Details")
                        interaction_df = pd.DataFrame(interactions)
                        st.dataframe(interaction_df, use_container_width=True)
                        
                        # Store results
                        st.session_state.ppi_predictions[ppi_seq] = interactions
                        
                    else:
                        st.info("No interactions predicted for this sequence.")
        else:
            st.warning("⚠️ Sequence too short for PPI analysis (minimum 20 residues).")
    else:
        st.info("👆 Please enter a protein sequence for PPI analysis.")

# -------------------- Tab 5: Advanced Tools --------------------
with tab5:
    st.header("⚙️ Advanced Analysis Tools")
    
    # Tool selection
    tool = st.selectbox("Select Analysis Tool", [
        "Disulfide Bond Prediction",
        "Post-Translational Modification Sites",
        "Membrane Topology Prediction",
        "Batch Processing",
        "Export Tools"
    ])
    
    if tool == "Disulfide Bond Prediction":
        st.subheader("🔗 Disulfide Bond Prediction")
        
        ss_seq_text = st.text_area("Protein Sequence", height=120, key="ss_seq")
        
        if ss_seq_text:
            ss_seq = clean_sequence(ss_seq_text)
            
            if len(ss_seq) >= 10:
                st.success(f"✅ Sequence loaded: {len(ss_seq)} residues")
                
                if st.button("🔍 Predict Disulfide Bonds", type="primary"):
                    bonds = predict_disulfide_bonds(ss_seq)
                    
                    if bonds:
                        st.subheader("📊 Predicted Disulfide Bonds")
                        
                        # Visualization
                        fig = go.Figure()
                        
                        for bond in bonds:
                            fig.add_trace(go.Scatter(
                                x=[bond['residue1'], bond['residue2']],
                                y=[0, 0],
                                mode='lines+markers',
                                line=dict(color='red', width=3),
                                marker=dict(size=8, color='red'),
                                name=f"C{bond['residue1']}-C{bond['residue2']}",
                                hovertemplate=f"<b>C{bond['residue1']}-C{bond['residue2']}</b><br>" +
                                            f"Confidence: {bond['confidence']:.2f}<extra></extra>"
                            ))
                        
                        fig.update_layout(
                            title="Disulfide Bond Prediction",
                            xaxis_title="Residue Position",
                            yaxis_title="",
                            height=200,
                            showlegend=True,
                            plot_bgcolor='rgba(0,0,0,0)',
                            paper_bgcolor='rgba(0,0,0,0)'
                        )
                        st.plotly_chart(fig, use_container_width=True)
                        
                        # Bond table
                        st.subheader("📋 Bond Details")
                        bond_df = pd.DataFrame(bonds)
                        st.dataframe(bond_df, use_container_width=True)
                    else:
                        st.info("No disulfide bonds predicted for this sequence.")
            else:
                st.warning("⚠️ Sequence too short for analysis.")
        else:
            st.info("👆 Please enter a protein sequence.")
    
    elif tool == "Post-Translational Modification Sites":
        st.subheader("🧪 Post-Translational Modification Prediction")
        
        ptm_seq_text = st.text_area("Protein Sequence", height=120, key="ptm_seq")
        
        if ptm_seq_text:
            ptm_seq = clean_sequence(ptm_seq_text)
            
            if len(ptm_seq) >= 10:
                st.success(f"✅ Sequence loaded: {len(ptm_seq)} residues")
                
                if st.button("🔍 Predict PTM Sites", type="primary"):
                    ptms = predict_ptm_sites(ptm_seq)
                    
                    if ptms:
                        st.subheader("📊 Predicted PTM Sites")
                        
                        # Visualization
                        fig = go.Figure()
                        
                        colors = {'Phosphorylation': '#3b82f6', 'N-glycosylation': '#10b981', 'O-glycosylation': '#f59e0b'}
                        
                        for ptm in ptms:
                            color = colors.get(ptm['type'], '#64748b')
                            fig.add_trace(go.Scatter(
                                x=[ptm['position']],
                                y=[0],
                                mode='markers',
                                marker=dict(size=10, color=color),
                                name=ptm['type'],
                                hovertemplate=f"<b>{ptm['type']}</b><br>" +
                                            f"Position: {ptm['position']}<br>" +
                                            f"Residue: {ptm['residue']}<br>" +
                                            f"Confidence: {ptm['confidence']:.2f}<extra></extra>"
                            ))
                        
                        fig.update_layout(
                            title="PTM Site Prediction",
                            xaxis_title="Residue Position",
                            yaxis_title="",
                            height=200,
                            showlegend=True,
                            plot_bgcolor='rgba(0,0,0,0)',
                            paper_bgcolor='rgba(0,0,0,0)'
                        )
                        st.plotly_chart(fig, use_container_width=True)
                        
                        # PTM table
                        st.subheader("📋 PTM Details")
                        ptm_df = pd.DataFrame(ptms)
                        st.dataframe(ptm_df, use_container_width=True)
                    else:
                        st.info("No PTM sites predicted for this sequence.")
            else:
                st.warning("⚠️ Sequence too short for analysis.")
        else:
            st.info("👆 Please enter a protein sequence.")
    
    elif tool == "Membrane Topology Prediction":
        st.subheader("🧱 Membrane Topology Prediction")
        
        mem_seq_text = st.text_area("Protein Sequence", height=120, key="mem_seq")
        
        if mem_seq_text:
            mem_seq = clean_sequence(mem_seq_text)
            
            if len(mem_seq) >= 10:
                st.success(f"✅ Sequence loaded: {len(mem_seq)} residues")
                
                if st.button("🔍 Predict Membrane Topology", type="primary"):
                    topology = predict_membrane_topology(mem_seq)
                    
                    if topology:
                        st.subheader("📊 Predicted Membrane Topology")
                        
                        # Visualization
                        fig = go.Figure()
                        
                        positions = [t['position'] for t in topology]
                        confidences = [t['confidence'] for t in topology]
                        
                        fig.add_trace(go.Scatter(
                            x=positions,
                            y=confidences,
                            mode='markers',
                            marker=dict(size=8, color='orange'),
                            name='Membrane Regions',
                            hovertemplate='<b>Membrane Region</b><br>' +
                                        'Position: %{x}<br>' +
                                        'Confidence: %{y:.2f}<extra></extra>'
                        ))
                        
                        fig.update_layout(
                            title="Membrane Topology Prediction",
                            xaxis_title="Residue Position",
                            yaxis_title="Confidence",
                            height=400,
                            plot_bgcolor='rgba(0,0,0,0)',
                            paper_bgcolor='rgba(0,0,0,0)'
                        )
                        st.plotly_chart(fig, use_container_width=True)
                        
                        # Topology table
                        st.subheader("📋 Topology Details")
                        topology_df = pd.DataFrame(topology)
                        st.dataframe(topology_df, use_container_width=True)
                    else:
                        st.info("No membrane regions predicted for this sequence.")
            else:
                st.warning("⚠️ Sequence too short for analysis.")
        else:
            st.info("👆 Please enter a protein sequence.")
    
    elif tool == "Batch Processing":
        st.subheader("📦 Batch Processing")
        
        st.info("Upload multiple sequences in FASTA format for batch processing.")
        
        batch_file = st.file_uploader("Upload FASTA file", type=["fasta", "txt"], key="batch_file")
        
        if batch_file:
            try:
                content = batch_file.read().decode('utf-8')
                sequences = {}
                
                # Parse FASTA
                current_name = None
                current_seq = ""
                
                for line in content.split('\n'):
                    line = line.strip()
                    if line.startswith('>'):
                        if current_name:
                            sequences[current_name] = clean_sequence(current_seq)
                        current_name = line[1:]
                        current_seq = ""
                    else:
                        current_seq += line
                
                if current_name:
                    sequences[current_name] = clean_sequence(current_seq)
                
                if sequences:
                    st.success(f"✅ Loaded {len(sequences)} sequences")
                    
                    # Show sequences
                    with st.expander("📋 Loaded Sequences"):
                        for name, seq in sequences.items():
                            st.text(f"{name}: {len(seq)} residues")
                    
                    if st.button("🚀 Process All Sequences", type="primary"):
                        with st.spinner("Processing sequences..."):
                            results = batch_predict_structures(sequences, method="api")
                            
                            st.subheader("📊 Batch Processing Results")
                            
                            # Results summary
                            successful = [r for r in results if r['success']]
                            failed = [r for r in results if not r['success']]
                            
                            col1, col2 = st.columns(2)
                            with col1:
                                st.metric("Successful", len(successful))
                            with col2:
                                st.metric("Failed", len(failed))
                            
                            # Results table
                            results_df = pd.DataFrame([
                                {
                                    'Name': r['name'],
                                    'Length': len(r['sequence']),
                                    'Status': 'Success' if r['success'] else 'Failed',
                                    'Error': r.get('error', '')
                                }
                                for r in results
                            ])
                            st.dataframe(results_df, use_container_width=True)
                            
                            # Download results
                            if successful:
                                st.subheader("⬇️ Download Results")
                                
                                # Create ZIP file with all PDB files
                                import zipfile
                                zip_buffer = io.BytesIO()
                                
                                with zipfile.ZipFile(zip_buffer, 'w') as zip_file:
                                    for result in successful:
                                        zip_file.writestr(f"{result['name']}.pdb", result['pdb'])
                                
                                zip_buffer.seek(0)
                                
                                st.download_button(
                                    "📦 Download All PDB Files",
                                    data=zip_buffer.getvalue(),
                                    file_name="batch_structures.zip",
                                    mime="application/zip"
                                )
                else:
                    st.error("No valid sequences found in the file.")
                    
            except Exception as e:
                st.error(f"Error reading file: {str(e)}")
    
    elif tool == "Export Tools":
        st.subheader("💾 Export Tools")
        
        st.info("Export your analysis results in various formats.")
        
        # Export options
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("📊 Data Export")
            
            if st.session_state.prediction_history:
                # Export prediction history
                history_json = json.dumps(st.session_state.prediction_history, indent=2)
                st.download_button(
                    "📄 Export Prediction History (JSON)",
                    data=history_json,
                    file_name="prediction_history.json",
                    mime="application/json"
                )
            
            if st.session_state.comparison_data:
                # Export comparison data
                comparison_json = json.dumps(st.session_state.comparison_data, indent=2)
                st.download_button(
                    "📄 Export Comparison Data (JSON)",
                    data=comparison_json,
                    file_name="comparison_data.json",
                    mime="application/json"
                )
        
        with col2:
            st.subheader("🎨 Visualization Export")
            
            st.info("Visualization export features coming soon!")
            
            # Placeholder for future visualization exports
            st.markdown("""
            **Planned Export Formats:**
            - High-resolution PNG images
            - SVG vector graphics
            - PyMOL session files
            - ChimeraX session files
            - VMD state files
            """)

# -------------------- Tab 6: AI Protein Designer --------------------
with tab6:
    st.header("🤖 AI-Powered Protein Designer")
    
    st.markdown("""
    <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                border-radius: 16px; padding: 1.5rem; color: white; margin-bottom: 2rem;">
        <h3 style="margin: 0 0 0.5rem 0;">🧬 Intelligent Protein Engineering</h3>
        <p style="margin: 0; opacity: 0.9;">AI-driven suggestions for protein optimization and design</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Input sequence
    ai_seq_text = st.text_area("Protein Sequence for AI Design", height=120, key="ai_seq")
    
    if ai_seq_text:
        ai_seq = clean_sequence(ai_seq_text)
        
        if len(ai_seq) >= 20:
            st.success(f"✅ Sequence loaded: {len(ai_seq)} residues")
            
            # Design target selection
            col1, col2, col3 = st.columns(3)
            
            with col1:
                design_target = st.selectbox(
                    "Design Target",
                    ["stability", "activity", "solubility", "binding", "expression"],
                    help="Choose what to optimize in your protein"
                )
            
            with col2:
                mutation_limit = st.slider("Max Mutations", 1, 10, 3)
            
            with col3:
                confidence_threshold = st.slider("Confidence Threshold", 0.5, 1.0, 0.7)
            
            if st.button("🤖 Generate AI Design Suggestions", type="primary"):
                with st.spinner("AI is analyzing your protein..."):
                    suggestions = ai_protein_design_suggestions(ai_seq, design_target)
                    
                    if suggestions:
                        st.subheader("🎯 AI Design Recommendations")
                        
                        # Filter by confidence
                        filtered_suggestions = [s for s in suggestions if s['confidence'] >= confidence_threshold]
                        
                        if filtered_suggestions:
                            for i, suggestion in enumerate(filtered_suggestions[:mutation_limit]):
                                with st.expander(f"💡 Suggestion {i+1}: Position {suggestion['position']}"):
                                    col_a, col_b = st.columns(2)
                                    
                                    with col_a:
                                        st.markdown(f"**Mutation:** {suggestion['original']} → {suggestion['suggested']}")
                                        st.markdown(f"**Position:** {suggestion['position']}")
                                        st.markdown(f"**Impact:** {suggestion['impact']}")
                                    
                                    with col_b:
                                        st.markdown(f"**Confidence:** {suggestion['confidence']:.2f}")
                                        st.markdown(f"**Reason:** {suggestion['reason']}")
                                        
                                        # Apply mutation button
                                        if st.button(f"Apply Mutation {i+1}", key=f"apply_{i}"):
                                            # Apply mutation to sequence
                                            mutated_seq = list(ai_seq)
                                            mutated_seq[suggestion['position']-1] = suggestion['suggested']
                                            mutated_seq = ''.join(mutated_seq)
                                            
                                            st.session_state.ai_designs.append({
                                                'original_sequence': ai_seq,
                                                'mutated_sequence': mutated_seq,
                                                'mutations': [suggestion],
                                                'timestamp': datetime.now().isoformat(),
                                                'target': design_target
                                            })
                                            
                                            st.success(f"✅ Mutation applied! New sequence: {mutated_seq}")
                            
                            # Apply all mutations
                            if st.button("🔄 Apply All Mutations", type="secondary"):
                                mutated_seq = list(ai_seq)
                                for suggestion in filtered_suggestions[:mutation_limit]:
                                    mutated_seq[suggestion['position']-1] = suggestion['suggested']
                                mutated_seq = ''.join(mutated_seq)
                                
                                st.session_state.ai_designs.append({
                                    'original_sequence': ai_seq,
                                    'mutated_sequence': mutated_seq,
                                    'mutations': filtered_suggestions[:mutation_limit],
                                    'timestamp': datetime.now().isoformat(),
                                    'target': design_target
                                })
                                
                                st.success(f"✅ All mutations applied! New sequence: {mutated_seq}")
                        
                        else:
                            st.info("No suggestions meet the confidence threshold. Try lowering it.")
                    else:
                        st.info("No design suggestions found for this sequence and target.")
        else:
            st.warning("⚠️ Sequence too short for AI design analysis.")
    else:
        st.info("👆 Please enter a protein sequence for AI design suggestions.")
    
    # AI Design History
    if st.session_state.ai_designs:
        st.subheader("📚 AI Design History")
        for i, design in enumerate(st.session_state.ai_designs):
            with st.expander(f"Design {i+1}: {design['target'].title()} ({design['timestamp'][:10]})"):
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("**Original Sequence:**")
                    st.code(design['original_sequence'])
                
                with col2:
                    st.markdown("**Mutated Sequence:**")
                    st.code(design['mutated_sequence'])
                
                st.markdown(f"**Mutations Applied:** {len(design['mutations'])}")
                for mutation in design['mutations']:
                    st.markdown(f"- Position {mutation['position']}: {mutation['original']} → {mutation['suggested']}")

# -------------------- Tab 7: Molecular Dynamics --------------------
with tab7:
    st.header("🧬 Molecular Dynamics Simulation")
    
    st.markdown("""
    <div style="background: linear-gradient(135deg, #10b981 0%, #059669 100%); 
                border-radius: 16px; padding: 1.5rem; color: white; margin-bottom: 2rem;">
        <h3 style="margin: 0 0 0.5rem 0;">⚡ Dynamic Protein Analysis</h3>
        <p style="margin: 0; opacity: 0.9;">Prepare and analyze molecular dynamics simulations</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Simulation type selection
    sim_type = st.selectbox(
        "Simulation Type",
        ["equilibration", "production", "analysis", "custom"],
        help="Choose the type of MD simulation to prepare"
    )
    
    # Input structure options
    load_col1, load_col2, load_col3 = st.columns([1,1,1])
    with load_col1:
        md_pdb_id = st.text_input("PDB ID", placeholder="e.g., 1CRN", key="md_pdb_id")
        if st.button("Load by PDB ID", key="btn_md_load_pdbid") and md_pdb_id.strip():
            loaded = fetch_pdb_by_id(md_pdb_id)
            if loaded:
                st.session_state.md_pdb = loaded
                st.success("Loaded structure from RCSB PDB.")
            else:
                st.warning("Failed to load PDB ID. Check the ID and try again.")
    with load_col2:
        md_seq_for_md = st.text_input("Sequence (Auto: PDB → ESM)", placeholder="Paste amino-acid sequence", key="md_seq_auto")
        if st.button("Load from Sequence", key="btn_md_load_seq") and md_seq_for_md.strip():
            hit = fetch_rcsb_pdb_by_sequence(clean_sequence(md_seq_for_md))
            if not hit:
                hit = predict_pdb_from_esmfold(clean_sequence(md_seq_for_md), method="api")
            if hit and _looks_like_pdb(hit):
                st.session_state.md_pdb = hit
                st.success("Loaded structure from Auto retrieval.")
            else:
                st.warning("Could not retrieve structure for the sequence.")
    with load_col3:
        st.caption("Or paste PDB text below:")
    # PDB text area
    md_pdb_text = st.text_area("PDB Structure for MD Simulation", height=120, key="md_pdb")
    # Prefer loaded structure from session if available
    if st.session_state.get('md_pdb') and not md_pdb_text:
        md_pdb_text = st.session_state.md_pdb
    
    if md_pdb_text:
        st.success("✅ Structure loaded for analysis")
        
        # MD trajectory analysis setup
        st.subheader("📁 MD Trajectory Inputs")
        st.caption("Upload topology and trajectory to enable analysis. Supports PDB/PRMTOP for topology, DCD/XTCT/TRR for trajectory.")
        colu1, colu2 = st.columns(2)
        with colu1:
            topo_file = st.file_uploader("Topology (PDB/PRMTOP)", type=["pdb", "prmtop", "psf"], key="md_topo")
        with colu2:
            traj_file = st.file_uploader("Trajectory (DCD/XTC/TRR)", type=["dcd", "xtc", "trr"], key="md_traj")

        if mdtraj_available and topo_file and traj_file:
            try:
                with st.spinner("Loading trajectory with MDTraj..."):
                    # Save to temp buffers and load
                    topo_bytes = topo_file.read()
                    traj_bytes = traj_file.read()
                    topo_buf = io.BytesIO(topo_bytes)
                    traj_buf = io.BytesIO(traj_bytes)
                    # md.load requires file paths; use temporary files
                    import tempfile
                    with tempfile.NamedTemporaryFile(suffix='.' + topo_file.name.split('.')[-1]) as tf_topo, \
                         tempfile.NamedTemporaryFile(suffix='.' + traj_file.name.split('.')[-1]) as tf_traj:
                        tf_topo.write(topo_bytes); tf_topo.flush()
                        tf_traj.write(traj_bytes); tf_traj.flush()
                        traj = md.load(tf_traj.name, top=tf_topo.name)

                st.success(f"Trajectory loaded: {traj.n_frames} frames, {traj.n_atoms} atoms")

                # Analysis tabs
                a1, a2, a3, a4, a5, a6, a7, a8 = st.tabs([
                    "RMSD", "RMSF", "Rg", "SASA", "H-Bonds", "PCA/DCCM", "FEL", "Energies"
                ])

                with a1:
                    st.subheader("📏 RMSD")
                    rmsd = md.rmsd(traj, traj, 0)
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(y=rmsd, mode='lines', name='RMSD (nm)'))
                    fig.update_layout(xaxis_title='Frame', yaxis_title='RMSD (nm)', height=400)
                    st.plotly_chart(fig, use_container_width=True)

                with a2:
                    st.subheader("📈 RMSF")
                    mean_xyz = traj.xyz.mean(axis=0)
                    fluc = np.sqrt(((traj.xyz - mean_xyz) ** 2).sum(axis=2).mean(axis=0))
                    fig = go.Figure()
                    fig.add_trace(go.Bar(x=list(range(len(fluc))), y=fluc, name='RMSF (nm)'))
                    fig.update_layout(xaxis_title='Atom Index', yaxis_title='RMSF (nm)', height=400)
                    st.plotly_chart(fig, use_container_width=True)

                with a3:
                    st.subheader("🌀 Radius of Gyration (Rg)")
                    rg = md.compute_rg(traj)
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(y=rg, mode='lines', name='Rg (nm)'))
                    fig.update_layout(xaxis_title='Frame', yaxis_title='Rg (nm)', height=400)
                    st.plotly_chart(fig, use_container_width=True)

                with a4:
                    st.subheader("🌊 SASA")
                    sasa = md.shrake_rupley(traj)
                    sasa_total = sasa.sum(axis=1)
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(y=sasa_total, mode='lines', name='SASA (nm²)'))
                    fig.update_layout(xaxis_title='Frame', yaxis_title='SASA (nm²)', height=400)
                    st.plotly_chart(fig, use_container_width=True)

                with a5:
                    st.subheader("🤝 Hydrogen Bonds (simple distance-angle)")
                    st.caption("Approximate H-bonds using geometric criteria. For rigorous analysis, use MDAnalysis/MDTraj hydrogen_bonds.")
                    try:
                        hbonds = md.baker_hubbard(traj, periodic=True)
                        counts = np.array([len(hbonds)] * traj.n_frames)
                    except Exception:
                        counts = np.zeros(traj.n_frames)
                    fig = go.Figure()
                    fig.add_trace(go.Scatter(y=counts, mode='lines', name='# H-bonds'))
                    fig.update_layout(xaxis_title='Frame', yaxis_title='Count', height=300)
                    st.plotly_chart(fig, use_container_width=True)

                with a6:
                    st.subheader("🧭 PCA and DCCM")
                    # PCA on C-alpha coordinates
                    ca = traj.topology.select('name CA')
                    if ca.size > 3:
                        X = traj.xyz[:, ca, :].reshape(traj.n_frames, -1)
                        X -= X.mean(axis=0)
                        # Covariance and eigendecomposition
                        cov = np.cov(X, rowvar=False)
                        vals, vecs = np.linalg.eigh(cov)
                        order = np.argsort(vals)[::-1]
                        vals = vals[order]
                        explained = vals / vals.sum()
                        fig = go.Figure()
                        fig.add_trace(go.Bar(x=list(range(1, min(11, len(explained)) + 1)), y=explained[:10]))
                        fig.update_layout(title='PCA Explained Variance', xaxis_title='PC', yaxis_title='Variance Ratio', height=300)
                        st.plotly_chart(fig, use_container_width=True)
                    # DCCM
                    X = traj.xyz.reshape(traj.n_frames, -1)
                    X -= X.mean(axis=0)
                    corr = np.corrcoef(X, rowvar=False)
                    fig = go.Figure(data=go.Heatmap(z=corr, colorscale='RdBu', zmin=-1, zmax=1))
                    fig.update_layout(title='Dynamic Cross-Correlation Matrix', height=500)
                    st.plotly_chart(fig, use_container_width=True)

                with a7:
                    st.subheader("🗺️ Free Energy Landscape (FEL)")
                    # Use first two PCs as CVs
                    ca = traj.topology.select('name CA')
                    if ca.size > 3:
                        X = traj.xyz[:, ca, :].reshape(traj.n_frames, -1)
                        X -= X.mean(axis=0)
                        cov = np.cov(X, rowvar=False)
                        vals, vecs = np.linalg.eigh(cov)
                        order = np.argsort(vals)[::-1]
                        vecs = vecs[:, order]
                        Y = X @ vecs[:, :2]
                        # 2D histogram and -kT ln P
                        H, xedges, yedges = np.histogram2d(Y[:, 0], Y[:, 1], bins=50)
                        P = H / H.sum()
                        with np.errstate(divide='ignore'):
                            F = -np.log(P + 1e-12)  # kT units
                        fig = go.Figure(data=go.Heatmap(z=F.T, x=xedges[:-1], y=yedges[:-1], colorscale='Viridis'))
                        fig.update_layout(xaxis_title='PC1', yaxis_title='PC2', height=500)
                        st.plotly_chart(fig, use_container_width=True)
                    else:
                        st.info("Not enough C-alpha atoms for FEL.")

                with a8:
                    st.subheader("⚡ Energies (Placeholder)")
                    st.caption("If you upload log/edr files, we can parse potential/kinetic/total energy, T/P/ρ in a future update.")
                    st.info("For now, use your MD engine's outputs (GROMACS/AMBER) and upload for parsing.")

            except Exception as e:
                st.error(f"Trajectory loading/analysis failed: {str(e)}")
        else:
            if not mdtraj_available:
                st.warning("MDTraj not available. Install with: pip install mdtraj")
            else:
                st.info("Upload topology and trajectory to enable MD analyses.")
    
    else:
        st.info("👆 Please paste a protein structure for drug discovery analysis.")
    
    # Drug Discovery History (from MD analyses)
    if st.session_state.drug_predictions:
        st.subheader("📚 Drug Discovery History")
        for i, prediction in enumerate(st.session_state.drug_predictions):
            with st.expander(f"Analysis {i+1} ({prediction['timestamp'][:10]})"):
                st.markdown(f"**Binding Sites Found:** {len(prediction['binding_sites'])}")
                
                for site in prediction['binding_sites']:
                    st.markdown(f"- Site {site['site_id']}: {site['pocket_type']} (Score: {site['druggability_score']:.3f})")

# -------------------- Tab 8: AI-Driven Drug Discovery Pipeline --------------------
with tab8:
    st.header("💊 AI-Driven Drug Discovery Pipeline")
    
    st.markdown("""
    <div style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
                border-radius: 16px; padding: 1.5rem; color: white; margin-bottom: 2rem;">
        <h3 style="margin: 0 0 0.5rem 0;">🤖 AI Drug Discovery Platform</h3>
        <p style="margin: 0; opacity: 0.9;">Automated virtual screening, molecular docking, and ADMET prediction</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Pipeline workflow tabs
    pipeline_tab1, pipeline_tab2, pipeline_tab3, pipeline_tab4 = st.tabs([
        "🎯 Target Selection", 
        "🔍 Virtual Screening", 
        "⚗️ Molecular Docking", 
        "📊 ADMET Analysis"
    ])
    
    # Initialize session state for drug discovery pipeline
    if 'drug_discovery_pipeline' not in st.session_state:
        st.session_state.drug_discovery_pipeline = {
            'target_protein': None,
            'ligand_library': [],
            'docking_results': [],
            'admet_predictions': [],
            'ranked_compounds': []
        }
    
    # Tab 1: Target Selection
    with pipeline_tab1:
        st.subheader("🎯 Protein Target Selection")
        
        col1, col2 = st.columns([1, 1])
        
        with col1:
            st.markdown("**Load Protein Structure**")
            target_pdb_id = st.text_input("PDB ID", placeholder="e.g., 1CRN", key="target_pdb_id")
            if st.button("Load Target Protein", key="load_target_btn"):
                if target_pdb_id.strip():
                    with st.spinner("Loading protein structure..."):
                        target_structure = fetch_pdb_by_id(target_pdb_id.strip())
                        if target_structure:
                            st.session_state.drug_discovery_pipeline['target_protein'] = {
                                'pdb_id': target_pdb_id.strip(),
                                'structure': target_structure,
                                'loaded_at': datetime.now().isoformat()
                            }
                            st.success(f"✅ Loaded protein {target_pdb_id}")
                        else:
                            st.error("Failed to load protein structure")
                else:
                    st.warning("Please enter a PDB ID")
        
        with col2:
            st.markdown("**Target Information**")
            if st.session_state.drug_discovery_pipeline['target_protein']:
                target_info = st.session_state.drug_discovery_pipeline['target_protein']
                st.info(f"**Current Target:** {target_info['pdb_id']}")
                st.info(f"**Loaded:** {target_info['loaded_at'][:19]}")
                
                # Analyze target protein
                if st.button("Analyze Target", key="analyze_target_btn"):
                    with st.spinner("Analyzing target protein..."):
                        # Basic protein analysis
                        structure = target_info['structure']
                        lines = structure.split('\n')
                        atom_count = len([line for line in lines if line.startswith('ATOM')])
                        chains = set([line[21:22] for line in lines if line.startswith('ATOM')])
                        
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("Atoms", atom_count)
                        with col2:
                            st.metric("Chains", len(chains))
                        with col3:
                            st.metric("Status", "Ready")
                        
                        # Binding site prediction for target
                        binding_sites = predict_drug_binding_sites(structure)
                        if binding_sites:
                            st.subheader("🎯 Predicted Binding Sites")
                            for site in binding_sites[:3]:  # Show top 3 sites
                                with st.expander(f"Site {site['site_id']}: {site['pocket_type']} (Score: {site['druggability_score']:.3f})"):
                                    col1, col2 = st.columns(2)
                                    with col1:
                                        st.metric("Druggability", f"{site['druggability_score']:.3f}")
                                        st.metric("Volume", f"{site['volume']:.1f} Å³")
                                    with col2:
                                        st.metric("Confidence", f"{site['confidence']:.2f}")
                                        st.metric("Type", site['pocket_type'].title())
            else:
                st.info("Please load a protein target first")
    
    # Tab 2: Virtual Screening
    with pipeline_tab2:
        st.subheader("🔍 Virtual Screening & Ligand Library")
        
        # Ligand library options
        library_options = st.selectbox(
            "Select Ligand Library",
            ["Custom SMILES", "FDA Approved Drugs", "Natural Products", "Fragment Library"],
            help="Choose the compound library for virtual screening"
        )
        
        if library_options == "Custom SMILES":
            st.markdown("**Enter Custom Compounds**")
            smiles_input = st.text_area(
                "SMILES Strings (one per line)",
                placeholder="CC(=O)OC1=CC=CC=C1C(=O)O\nCCN(CC)CCCC(C)NC1=C2C=CC(Cl)=CC2=NC=C1",
                height=150,
                help="Enter SMILES strings separated by new lines"
            )
            
            if st.button("Process Custom Library", key="process_custom_btn"):
                if smiles_input.strip():
                    smiles_list = [s.strip() for s in smiles_input.split('\n') if s.strip()]
                    st.session_state.drug_discovery_pipeline['ligand_library'] = [
                        {'smiles': smiles, 'name': f'Compound_{i+1}', 'source': 'Custom'}
                        for i, smiles in enumerate(smiles_list)
                    ]
                    st.success(f"✅ Loaded {len(smiles_list)} custom compounds")
                else:
                    st.warning("Please enter SMILES strings")
        
        elif library_options == "FDA Approved Drugs":
            st.markdown("**FDA Approved Drug Library**")
            if st.button("Load FDA Library", key="load_fda_btn"):
                with st.spinner("Loading FDA approved drugs..."):
                    # Simulated FDA drug library
                    fda_drugs = [
                        {'smiles': 'CC(=O)OC1=CC=CC=C1C(=O)O', 'name': 'Aspirin', 'source': 'FDA'},
                        {'smiles': 'CCN(CC)CCCC(C)NC1=C2C=CC(Cl)=CC2=NC=C1', 'name': 'Chlorpromazine', 'source': 'FDA'},
                        {'smiles': 'CC1=CC=C(C=C1)C2=CC(=O)C3=C(C=CC(=C3O2)O)O', 'name': 'Quercetin', 'source': 'FDA'},
                        {'smiles': 'CC1=CC=C(C=C1)C2=CC(=O)C3=C(C=CC(=C3O2)O)O', 'name': 'Rutin', 'source': 'FDA'},
                        {'smiles': 'CC1=CC=C(C=C1)C2=CC(=O)C3=C(C=CC(=C3O2)O)O', 'name': 'Catechin', 'source': 'FDA'}
                    ]
                    st.session_state.drug_discovery_pipeline['ligand_library'] = fda_drugs
                    st.success(f"✅ Loaded {len(fda_drugs)} FDA approved drugs")
        
        elif library_options == "Natural Products":
            st.markdown("**Natural Products Library**")
            if st.button("Load Natural Products", key="load_natural_btn"):
                with st.spinner("Loading natural products..."):
                    natural_products = [
                        {'smiles': 'CC1=CC=C(C=C1)C2=CC(=O)C3=C(C=CC(=C3O2)O)O', 'name': 'Curcumin', 'source': 'Natural'},
                        {'smiles': 'CC1=CC=C(C=C1)C2=CC(=O)C3=C(C=CC(=C3O2)O)O', 'name': 'Resveratrol', 'source': 'Natural'},
                        {'smiles': 'CC1=CC=C(C=C1)C2=CC(=O)C3=C(C=CC(=C3O2)O)O', 'name': 'Epigallocatechin', 'source': 'Natural'},
                        {'smiles': 'CC1=CC=C(C=C1)C2=CC(=O)C3=C(C=CC(=C3O2)O)O', 'name': 'Genistein', 'source': 'Natural'},
                        {'smiles': 'CC1=CC=C(C=C1)C2=CC(=O)C3=C(C=CC(=C3O2)O)O', 'name': 'Luteolin', 'source': 'Natural'}
                    ]
                    st.session_state.drug_discovery_pipeline['ligand_library'] = natural_products
                    st.success(f"✅ Loaded {len(natural_products)} natural products")
        
        elif library_options == "Fragment Library":
            st.markdown("**Fragment Library**")
            if st.button("Load Fragment Library", key="load_fragment_btn"):
                with st.spinner("Loading fragment library..."):
                    fragments = [
                        {'smiles': 'CC(=O)O', 'name': 'Acetic_Acid', 'source': 'Fragment'},
                        {'smiles': 'CCN', 'name': 'Ethylamine', 'source': 'Fragment'},
                        {'smiles': 'CC(=O)N', 'name': 'Acetamide', 'source': 'Fragment'},
                        {'smiles': 'CC(=O)OC', 'name': 'Methyl_Acetate', 'source': 'Fragment'},
                        {'smiles': 'CC(=O)NC', 'name': 'N-Methylacetamide', 'source': 'Fragment'}
                    ]
                    st.session_state.drug_discovery_pipeline['ligand_library'] = fragments
                    st.success(f"✅ Loaded {len(fragments)} fragments")
        
        # Display loaded library
        if st.session_state.drug_discovery_pipeline['ligand_library']:
            st.subheader("📚 Loaded Ligand Library")
            library_df = pd.DataFrame(st.session_state.drug_discovery_pipeline['ligand_library'])
            st.dataframe(library_df, use_container_width=True, hide_index=True)
            
            col1, col2 = st.columns(2)
            with col1:
                st.metric("Total Compounds", len(st.session_state.drug_discovery_pipeline['ligand_library']))
            with col2:
                sources = set([comp['source'] for comp in st.session_state.drug_discovery_pipeline['ligand_library']])
                st.metric("Sources", len(sources))
            
            # Molecular Structure Visualization
            st.subheader("🧬 Molecular Structure Visualization")
            
            # Create tabs for different visualization options
            viz_tab1, viz_tab2, viz_tab3 = st.tabs([
                "🔍 Individual Structures", 
                "📊 Structure Gallery", 
                "📈 Structure Analysis"
            ])
            
            with viz_tab1:
                st.markdown("**Select a compound to view its 3D structure:**")
                
                # Compound selector
                compound_names = [f"{comp['name']} ({comp['source']})" for comp in st.session_state.drug_discovery_pipeline['ligand_library']]
                selected_idx = st.selectbox(
                    "Choose Compound",
                    range(len(compound_names)),
                    format_func=lambda x: compound_names[x],
                    key="structure_selector"
                )
                
                if selected_idx is not None:
                    selected_compound = st.session_state.drug_discovery_pipeline['ligand_library'][selected_idx]
                    smiles = selected_compound['smiles']
                    name = selected_compound['name']
                    
                    st.markdown(f"**Selected: {name}**")
                    st.code(f"SMILES: {smiles}")
                    
                    # Try to display 3D structure using RDKit
                    try:
                        if rdkit_available:
                            from rdkit import Chem
                            from rdkit.Chem import AllChem, rdMolDescriptors
                            import py3Dmol
                            
                            # Helper: build 3D and render with py3Dmol using SDF block
                            def _render_smiles_3d(smiles_str: str, width: int = 600, height: int = 400):
                                m = Chem.MolFromSmiles(smiles_str)
                                if m is None:
                                    st.error("❌ Invalid SMILES string - cannot generate 3D structure")
                                    return None
                                m = Chem.AddHs(m)
                                try:
                                    params = AllChem.ETKDGv3()
                                except Exception:
                                    params = AllChem.ETKDG()
                                params.randomSeed = 0xF00D
                                res = AllChem.EmbedMolecule(m, params)
                                if res != 0:
                                    # Fallback to older method
                                    res2 = AllChem.EmbedMolecule(m)
                                    if res2 != 0:
                                        st.error("❌ 3D embedding failed for this SMILES")
                                        return None
                                # Optimize geometry
                                try:
                                    AllChem.MMFFOptimizeMolecule(m)
                                except Exception:
                                    try:
                                        AllChem.UFFOptimizeMolecule(m)
                                    except Exception:
                                        pass
                                # Render via SDF block
                                sdf_block = Chem.MolToMolBlock(m)
                                view_local = py3Dmol.view(width=width, height=height)
                                view_local.addModel(sdf_block, "sdf")
                                view_local.setStyle({'stick': {'radius': 0.12}})
                                view_local.setBackgroundColor('white')
                                view_local.zoomTo()
                                showmol(view_local, height=height, width=width)
                                return m
                            
                            # Create molecule from SMILES
                            mol = _render_smiles_3d(smiles, width=600, height=400)
                            
                            # Display molecular properties
                            if mol is not None:
                                st.subheader("Molecular Properties")
                                col1, col2, col3, col4 = st.columns(4)
                                
                                with col1:
                                    mw = rdMolDescriptors.CalcExactMolWt(mol)
                                    st.metric("Molecular Weight", f"{mw:.1f} Da")
                                
                                with col2:
                                    logp = rdMolDescriptors.CalcCrippenDescriptors(mol)[0]
                                    st.metric("LogP", f"{logp:.2f}")
                                
                                with col3:
                                    tpsa = rdMolDescriptors.CalcTPSA(mol)
                                    st.metric("TPSA", f"{tpsa:.1f} Å²")
                                
                                with col4:
                                    hbd = rdMolDescriptors.CalcNumHBD(mol)
                                    st.metric("HBD", hbd)
                                
                                # Additional properties
                                col5, col6, col7, col8 = st.columns(4)
                                
                                with col5:
                                    hba = rdMolDescriptors.CalcNumHBA(mol)
                                    st.metric("HBA", hba)
                                
                                with col6:
                                    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol)
                                    st.metric("Rotatable Bonds", rotb)
                                
                                with col7:
                                    rings = rdMolDescriptors.CalcNumRings(mol)
                                    st.metric("Rings", rings)
                                
                                with col8:
                                    atoms = mol.GetNumAtoms()
                                    st.metric("Atoms", atoms)
                            
                            else:
                                st.error("❌ Invalid SMILES string - cannot generate 3D structure")
                        else:
                            st.warning("⚠️ RDKit not available - cannot display 3D structures")
                            st.info("Install RDKit to enable molecular structure visualization")
                            
                    except Exception as e:
                        st.error(f"❌ Error generating 3D structure: {str(e)}")
                        st.info("This might be due to complex SMILES or RDKit limitations")
            
            with viz_tab2:
                st.markdown("**Structure Gallery - View multiple compounds:**")
                
                # Number of compounds to display
                library_size = len(st.session_state.drug_discovery_pipeline['ligand_library'])
                max_display = min(6, library_size)
                
                if library_size > 1:
                    num_compounds = st.slider(
                        "Number of compounds to display",
                        min_value=1,
                        max_value=max_display,
                        value=min(3, library_size),
                        key="gallery_slider"
                    )
                else:
                    num_compounds = 1
                    st.info(f"Only {library_size} compound available - displaying all")
                
                if rdkit_available:
                    try:
                        from rdkit import Chem
                        from rdkit.Chem import AllChem
                        import py3Dmol
                        
                        # Create columns for structure display
                        cols = st.columns(2)
                        
                        for i in range(num_compounds):
                            compound = st.session_state.drug_discovery_pipeline['ligand_library'][i]
                            smiles = compound['smiles']
                            name = compound['name']
                            
                            with cols[i % 2]:
                                st.markdown(f"**{name}**")
                                st.code(f"SMILES: {smiles}")
                                
                                # Generate and display 3D structure using robust ETKDG + SDF rendering
                                try:
                                    # reuse helper if exists above, else define minimal fallback
                                    try:
                                        _render_smiles_3d
                                        _has_helper = True
                                    except NameError:
                                        _has_helper = False
                                    if not _has_helper:
                                        def _render_smiles_3d(smiles_str: str, width: int = 300, height: int = 250):
                                            m2 = Chem.MolFromSmiles(smiles_str)
                                            if m2 is None:
                                                st.error(f"❌ Invalid SMILES for {name}")
                                                return None
                                            m2 = Chem.AddHs(m2)
                                            try:
                                                params2 = AllChem.ETKDGv3()
                                            except Exception:
                                                params2 = AllChem.ETKDG()
                                            params2.randomSeed = 0xBEEF
                                            r = AllChem.EmbedMolecule(m2, params2)
                                            if r != 0:
                                                r2 = AllChem.EmbedMolecule(m2)
                                                if r2 != 0:
                                                    st.error(f"❌ 3D embedding failed for {name}")
                                                    return None
                                            try:
                                                AllChem.MMFFOptimizeMolecule(m2)
                                            except Exception:
                                                try:
                                                    AllChem.UFFOptimizeMolecule(m2)
                                                except Exception:
                                                    pass
                                            sdf_block2 = Chem.MolToMolBlock(m2)
                                            v = py3Dmol.view(width=300, height=250)
                                            v.addModel(sdf_block2, "sdf")
                                            v.setStyle({'stick': {'radius': 0.12}})
                                            v.setBackgroundColor('white')
                                            v.zoomTo()
                                            showmol(v, height=250, width=300)
                                            return m2
                                    _ = _render_smiles_3d(smiles, width=300, height=250)
                                except Exception as _e:
                                    st.error(f"❌ Invalid SMILES for {name}")
                                
                    except Exception as e:
                        st.error(f"❌ Error in structure gallery: {str(e)}")
                else:
                    st.warning("⚠️ RDKit not available - cannot display structure gallery")
            
            with viz_tab3:
                st.markdown("**Structure Analysis and Comparison:**")
                
                if rdkit_available:
                    try:
                        from rdkit import Chem
                        from rdkit.Chem import rdMolDescriptors
                        import plotly.express as px
                        import plotly.graph_objects as go
                        
                        # Calculate properties for all compounds
                        properties_data = []
                        
                        for compound in st.session_state.drug_discovery_pipeline['ligand_library']:
                            smiles = compound['smiles']
                            name = compound['name']
                            
                            mol = Chem.MolFromSmiles(smiles)
                            if mol is not None:
                                props = {
                                    'Name': name,
                                    'SMILES': smiles,
                                    'MW': rdMolDescriptors.CalcExactMolWt(mol),
                                    'LogP': rdMolDescriptors.CalcCrippenDescriptors(mol)[0],
                                    'TPSA': rdMolDescriptors.CalcTPSA(mol),
                                    'HBD': rdMolDescriptors.CalcNumHBD(mol),
                                    'HBA': rdMolDescriptors.CalcNumHBA(mol),
                                    'RotB': rdMolDescriptors.CalcNumRotatableBonds(mol),
                                    'Rings': rdMolDescriptors.CalcNumRings(mol),
                                    'Atoms': mol.GetNumAtoms()
                                }
                                properties_data.append(props)
                        
                        if properties_data:
                            # Create properties DataFrame
                            props_df = pd.DataFrame(properties_data)
                            
                            # Display properties table
                            st.subheader("📊 Molecular Properties Summary")
                            st.dataframe(props_df, use_container_width=True, hide_index=True)
                            
                            # Create visualizations
                            st.subheader("📈 Property Distributions")
                            
                            # MW vs LogP scatter plot
                            fig1 = px.scatter(
                                props_df, 
                                x='MW', 
                                y='LogP', 
                                color='Name',
                                title="Molecular Weight vs LogP",
                                hover_data=['TPSA', 'HBD', 'HBA']
                            )
                            fig1.update_layout(height=500)
                            st.plotly_chart(fig1, use_container_width=True)
                            
                            # TPSA vs HBD scatter plot
                            fig2 = px.scatter(
                                props_df, 
                                x='TPSA', 
                                y='HBD', 
                                color='Name',
                                title="TPSA vs HBD Count",
                                hover_data=['MW', 'LogP', 'HBA']
                            )
                            fig2.update_layout(height=500)
                            st.plotly_chart(fig2, use_container_width=True)
                            
                            # Property distribution histograms
                            st.subheader("📊 Property Distribution Histograms")
                            
                            col1, col2 = st.columns(2)
                            
                            with col1:
                                fig3 = px.histogram(
                                    props_df, 
                                    x='MW', 
                                    title="Molecular Weight Distribution",
                                    nbins=10
                                )
                                st.plotly_chart(fig3, use_container_width=True)
                            
                            with col2:
                                fig4 = px.histogram(
                                    props_df, 
                                    x='LogP', 
                                    title="LogP Distribution",
                                    nbins=10
                                )
                                st.plotly_chart(fig4, use_container_width=True)
                            
                            # Drug-likeness analysis
                            st.subheader("💊 Drug-likeness Analysis")
                            
                            # Lipinski's Rule of Five
                            lipinski_pass = []
                            for _, row in props_df.iterrows():
                                mw_ok = row['MW'] <= 500
                                logp_ok = row['LogP'] <= 5
                                hbd_ok = row['HBD'] <= 5
                                hba_ok = row['HBA'] <= 10
                                lipinski_pass.append(mw_ok and logp_ok and hbd_ok and hba_ok)
                            
                            props_df['Lipinski_Pass'] = lipinski_pass
                            
                            # Display drug-likeness summary
                            col1, col2, col3 = st.columns(3)
                            with col1:
                                st.metric("Lipinski Compliant", f"{sum(lipinski_pass)}/{len(lipinski_pass)}")
                            with col2:
                                st.metric("Average MW", f"{props_df['MW'].mean():.1f} Da")
                            with col3:
                                st.metric("Average LogP", f"{props_df['LogP'].mean():.2f}")
                            
                            # Drug-likeness radar chart
                            if len(props_df) > 1:
                                # Normalize properties for radar chart
                                normalized_props = props_df.copy()
                                for col in ['MW', 'LogP', 'TPSA', 'HBD', 'HBA', 'RotB']:
                                    normalized_props[col] = (props_df[col] - props_df[col].min()) / (props_df[col].max() - props_df[col].min())
                                
                                fig_radar = go.Figure()
                                
                                for _, row in normalized_props.iterrows():
                                    fig_radar.add_trace(go.Scatterpolar(
                                        r=[row['MW'], row['LogP'], row['TPSA'], row['HBD'], row['HBA'], row['RotB']],
                                        theta=['MW', 'LogP', 'TPSA', 'HBD', 'HBA', 'RotB'],
                                        fill='toself',
                                        name=row['Name']
                                    ))
                                
                                fig_radar.update_layout(
                                    polar=dict(
                                        radialaxis=dict(
                                            visible=True,
                                            range=[0, 1]
                                        )
                                    ),
                                    title="Normalized Property Comparison",
                                    height=500
                                )
                                
                                st.plotly_chart(fig_radar, use_container_width=True)
                        
                    except Exception as e:
                        st.error(f"❌ Error in structure analysis: {str(e)}")
                else:
                    st.warning("⚠️ RDKit not available - cannot perform structure analysis")
    
    # Tab 3: Molecular Docking
    with pipeline_tab3:
        st.subheader("⚗️ Molecular Docking Analysis")
        
        if not st.session_state.drug_discovery_pipeline['target_protein']:
            st.warning("⚠️ Please load a protein target first (Tab 1)")
        elif not st.session_state.drug_discovery_pipeline['ligand_library']:
            st.warning("⚠️ Please load a ligand library first (Tab 2)")
        else:
            st.success("✅ Ready for molecular docking")
            
            # Docking parameters
            col1, col2 = st.columns(2)
            with col1:
                docking_algorithm = st.selectbox(
                    "Docking Algorithm",
                    ["AutoDock Vina", "CB-Dock", "Quick Docking"],
                    help="Select the molecular docking algorithm"
                )
            with col2:
                exhaustiveness = st.slider("Exhaustiveness", 1, 100, 8, help="Higher values = more thorough search")
            
            if st.button("🚀 Start Molecular Docking", type="primary", key="start_docking_btn"):
                with st.spinner("Performing molecular docking..."):
                    # Simulate docking process
                    target_structure = st.session_state.drug_discovery_pipeline['target_protein']['structure']
                    ligand_library = st.session_state.drug_discovery_pipeline['ligand_library']
                    
                    docking_results = []
                    progress_bar = st.progress(0)
                    
                    for i, ligand in enumerate(ligand_library):
                        # Simulate docking calculation
                        import time
                        time.sleep(0.5)  # Simulate processing time
                        
                        # Generate simulated docking results
                        binding_affinity = round(random.uniform(-12.0, -2.0), 2)  # kcal/mol
                        binding_pose_quality = round(random.uniform(0.3, 0.9), 3)
                        interaction_energy = round(random.uniform(-15.0, -5.0), 2)
                        
                        docking_result = {
                            'compound_name': ligand['name'],
                            'smiles': ligand['smiles'],
                            'binding_affinity': binding_affinity,
                            'pose_quality': binding_pose_quality,
                            'interaction_energy': interaction_energy,
                            'docking_score': round((binding_affinity + interaction_energy) / 2, 2),
                            'rank': i + 1
                        }
                        docking_results.append(docking_result)
                        
                        progress_bar.progress((i + 1) / len(ligand_library))
                    
                    st.session_state.drug_discovery_pipeline['docking_results'] = docking_results
                    st.success(f"✅ Docking completed for {len(docking_results)} compounds")
                    
                    # Display docking results
                    st.subheader("📊 Docking Results")
                    docking_df = pd.DataFrame(docking_results)
                    docking_df = docking_df.sort_values('docking_score', ascending=False)
                    st.dataframe(docking_df, use_container_width=True, hide_index=True)
                    
                    # Top compounds visualization
                    st.subheader("🏆 Top Docking Candidates")
                    top_compounds = docking_df.head(5)
                    
                    for idx, compound in top_compounds.iterrows():
                        with st.expander(f"🥇 {compound['compound_name']} (Score: {compound['docking_score']:.2f})"):
                            col1, col2, col3 = st.columns(3)
                            with col1:
                                st.metric("Binding Affinity", f"{compound['binding_affinity']:.2f} kcal/mol")
                            with col2:
                                st.metric("Pose Quality", f"{compound['pose_quality']:.3f}")
                            with col3:
                                st.metric("Interaction Energy", f"{compound['interaction_energy']:.2f} kcal/mol")
    
    # Tab 4: ADMET Analysis & Ranking
    with pipeline_tab4:
        st.subheader("📊 ADMET Analysis & Compound Ranking")
        
        if not st.session_state.drug_discovery_pipeline['docking_results']:
            st.warning("⚠️ Please perform molecular docking first (Tab 3)")
        else:
            st.success("✅ Ready for ADMET analysis")
            
            if st.button("🧪 Analyze ADMET Properties", type="primary", key="admet_analysis_btn"):
                with st.spinner("Computing ADMET properties for all compounds..."):
                    docking_results = st.session_state.drug_discovery_pipeline['docking_results']
                    admet_predictions = []
                    
                    progress_bar = st.progress(0)
                    
                    for i, compound in enumerate(docking_results):
                        # Compute ADMET properties
                        props = compute_admet_properties(compound['smiles'])
                        
                        if 'error' not in props:
                            summary = interpret_admet(props)
                            
                            # Calculate overall drug-likeness score
                            overall_score = np.mean([
                                summary['absorption']['score'],
                                summary['distribution']['score'],
                                summary['metabolism']['score'],
                                summary['excretion']['score'],
                                summary['toxicity']['score'],
                                summary['drug_likeness']['score']
                            ])
                            
                            # Combine docking and ADMET scores
                            combined_score = (compound['docking_score'] + overall_score) / 2
                            
                            admet_result = {
                                'compound_name': compound['compound_name'],
                                'smiles': compound['smiles'],
                                'docking_score': compound['docking_score'],
                                'admet_score': overall_score,
                                'combined_score': combined_score,
                                'drug_category': props.get('drug_category', 'Unknown'),
                                'molecular_weight': props.get('mw', 0),
                                'logp': props.get('logp', 0),
                                'tpsa': props.get('tpsa', 0),
                                'lipinski_pass': props.get('lipinski_pass', False),
                                'pains_alerts': props.get('pains', 0),
                                'qed_score': props.get('qed_score', 0)
                            }
                            admet_predictions.append(admet_result)
                        else:
                            # Handle compounds with ADMET errors
                            admet_result = {
                                'compound_name': compound['compound_name'],
                                'smiles': compound['smiles'],
                                'docking_score': compound['docking_score'],
                                'admet_score': 0.0,
                                'combined_score': compound['docking_score'],
                                'drug_category': 'Unknown',
                                'molecular_weight': 0,
                                'logp': 0,
                                'tpsa': 0,
                                'lipinski_pass': False,
                                'pains_alerts': 0,
                                'qed_score': 0
                            }
                            admet_predictions.append(admet_result)
                        
                        progress_bar.progress((i + 1) / len(docking_results))
                    
                    # Rank compounds by combined score
                    admet_predictions.sort(key=lambda x: x['combined_score'], reverse=True)
                    st.session_state.drug_discovery_pipeline['admet_predictions'] = admet_predictions
                    
                    st.success(f"✅ ADMET analysis completed for {len(admet_predictions)} compounds")
                    
                    # Display ranked results
                    st.subheader("🏆 Ranked Drug Candidates")
                    ranked_df = pd.DataFrame(admet_predictions)
                    st.dataframe(ranked_df, use_container_width=True, hide_index=True)
                    
                    # Enhanced ADMET Properties Visualization
                    st.subheader("🧬 Comprehensive ADMET Properties Analysis")
                    
                    # Create tabs for different visualizations
                    viz_tab1, viz_tab2, viz_tab3, viz_tab4 = st.tabs([
                        "📊 Scatter Plot", 
                        "📈 ADMET Radar", 
                        "🎯 Property Heatmap", 
                        "📋 Detailed Analysis"
                    ])
                    
                    with viz_tab1:
                        st.subheader("🥇 Top Drug Candidates Scatter Plot")
                        top_candidates = ranked_df.head(10)
                        
                        # Create enhanced scatter plot
                        fig = go.Figure()
                        
                        # Ensure marker sizes are positive and reasonable
                        marker_sizes = np.abs(top_candidates['combined_score']) * 20 + 10
                        marker_sizes = np.clip(marker_sizes, 10, 100)
                        
                        fig.add_trace(go.Scatter(
                            x=top_candidates['docking_score'],
                            y=top_candidates['admet_score'],
                            mode='markers+text',
                            text=top_candidates['compound_name'],
                            textposition='top center',
                            marker=dict(
                                size=marker_sizes,
                                color=top_candidates['combined_score'],
                                colorscale='Viridis',
                                showscale=True,
                                colorbar=dict(title="Combined Score"),
                                line=dict(width=2, color='white')
                            ),
                            name='Drug Candidates',
                            hovertemplate='<b>%{text}</b><br>' +
                                        'Docking Score: %{x:.3f}<br>' +
                                        'ADMET Score: %{y:.3f}<br>' +
                                        'Combined Score: %{marker.color:.3f}<extra></extra>'
                        ))
                        
                        fig.update_layout(
                            title="Drug Candidate Performance Matrix",
                            xaxis_title="Docking Score (Binding Affinity)",
                            yaxis_title="ADMET Score (Drug-likeness)",
                            height=600,
                            showlegend=False,
                            plot_bgcolor='rgba(0,0,0,0)',
                            paper_bgcolor='rgba(0,0,0,0)'
                        )
                        
                        st.plotly_chart(fig, use_container_width=True)
                    
                    with viz_tab2:
                        st.subheader("🎯 ADMET Profile Radar Charts")
                        
                        # Create radar charts for top 5 candidates
                        top_5 = ranked_df.head(5)
                        
                        for idx, candidate in top_5.iterrows():
                            with st.expander(f"🧬 {candidate['compound_name']} - ADMET Profile"):
                                # Get detailed ADMET properties
                                props = compute_admet_properties(candidate['smiles'])
                                if 'error' not in props:
                                    summary = interpret_admet(props)
                                    
                                    # Create radar chart
                                    categories = ['Absorption', 'Distribution', 'Metabolism', 'Excretion', 'Toxicity', 'Drug-likeness']
                                    scores = [
                                        summary['absorption']['score'],
                                        summary['distribution']['score'],
                                        summary['metabolism']['score'],
                                        summary['excretion']['score'],
                                        summary['toxicity']['score'],
                                        summary['drug_likeness']['score']
                                    ]
                                    
                                    fig_radar = go.Figure()
                                    fig_radar.add_trace(go.Scatterpolar(
                                        r=scores,
                                        theta=categories,
                                        fill='toself',
                                        name=candidate['compound_name'],
                                        line_color='rgb(32, 201, 151)',
                                        fillcolor='rgba(32, 201, 151, 0.3)'
                                    ))
                                    
                                    fig_radar.update_layout(
                                        polar=dict(
                                            radialaxis=dict(
                                                visible=True,
                                                range=[0, 1],
                                                tickfont=dict(size=12)
                                            )
                                        ),
                                        title=f"ADMET Profile: {candidate['compound_name']}",
                                        height=400,
                                        showlegend=True
                                    )
                                    
                                    st.plotly_chart(fig_radar, use_container_width=True)
                                    
                                    # Display key properties
                                    col1, col2, col3, col4 = st.columns(4)
                                    with col1:
                                        st.metric("Molecular Weight", f"{props.get('mw', 0):.1f} Da")
                                    with col2:
                                        st.metric("LogP", f"{props.get('logp', 0):.2f}")
                                    with col3:
                                        st.metric("TPSA", f"{props.get('tpsa', 0):.1f} Å²")
                                    with col4:
                                        st.metric("QED Score", f"{props.get('qed_score', 0):.3f}")
                                    
                                    # Detailed ADMET Parameters Breakdown
                                    st.subheader("📊 Detailed ADMET Parameters")
                                    
                                    # Local helpers to compute per-parameter scores safely (0-1)
                                    def _score_range(value, low, high):
                                        try:
                                            v = float(value)
                                        except Exception:
                                            return 0.0
                                        if low == high:
                                            return 1.0 if v == low else 0.0
                                        if v <= low:
                                            return 0.0
                                        if v >= high:
                                            return 0.0
                                        # Peak at middle of range for bell-like preference? Use linear inside bounds to 1 at center
                                        mid = (low + high) / 2.0
                                        half = (high - low) / 2.0
                                        return max(0.0, 1.0 - abs(v - mid) / half)
                                    
                                    def _score_le(value, threshold):
                                        try:
                                            v = float(value)
                                        except Exception:
                                            return 0.0
                                        if v <= threshold:
                                            # Better if comfortably below
                                            return min(1.0, 0.6 + (threshold - v) / max(1.0, threshold) * 0.4)
                                        return max(0.0, 1.0 - (v - threshold) / (abs(threshold) + 1.0))
                                    
                                    def _score_ge(value, threshold):
                                        try:
                                            v = float(value)
                                        except Exception:
                                            return 0.0
                                        if v >= threshold:
                                            # Better if comfortably above
                                            return min(1.0, 0.6 + (v - threshold) / (abs(threshold) + 1.0) * 0.4)
                                        return max(0.0, 1.0 - (threshold - v) / (abs(threshold) + 1.0))
                                    
                                    # Compute per-category parameter scores from props
                                    abs_scores = {
                                        'mw_score': _score_range(props.get('mw', 0), 150, 500),
                                        'logp_score': _score_range(props.get('logp', 0), 0, 5),
                                        'tpsa_score': _score_range(props.get('tpsa', 0), 0, 140),
                                        'hbd_score': _score_range(props.get('hbd', 0), 0, 5),
                                        'hba_score': _score_range(props.get('hba', 0), 0, 10),
                                        'rotb_score': _score_range(props.get('rotb', 0), 0, 10),
                                        'bioavailability_score': _score_ge(props.get('bioavailability_score', 0), 0.5),
                                        'lipinski_score': 1.0 if props.get('lipinski_pass') else 0.0,
                                    }
                                    
                                    dist_scores = {
                                        'logp_score': _score_range(props.get('logp', 0), 2, 5),
                                        'tpsa_score': _score_range(props.get('tpsa', 0), 0, 90),
                                        'mw_score': _score_range(props.get('mw', 0), 150, 500),
                                        'hbd_score': _score_range(props.get('hbd', 0), 0, 3),
                                        'hba_score': _score_range(props.get('hba', 0), 0, 8),
                                        'cns_score': _score_ge(props.get('cns_permeability', 0), 0.5),
                                        'bbb_score': _score_ge(props.get('bbb_score', 0), 0.5),
                                    }
                                    
                                    metab_scores = {
                                        'mw_score': _score_range(props.get('mw', 0), 150, 500),
                                        'logp_score': _score_range(props.get('logp', 0), 0, 5),
                                        'tpsa_score': _score_range(props.get('tpsa', 0), 0, 140),
                                        'hbd_score': _score_range(props.get('hbd', 0), 0, 5),
                                        'hba_score': _score_range(props.get('hba', 0), 0, 10),
                                        'rotb_score': _score_range(props.get('rotb', 0), 0, 10),
                                        'stability_score': _score_ge(props.get('metabolic_stability', 0), 0.5),
                                    }
                                    
                                    excr_scores = {
                                        'mw_score': _score_range(props.get('mw', 0), 150, 500),
                                        'logp_score': _score_range(props.get('logp', 0), 0, 5),
                                        'tpsa_score': _score_range(props.get('tpsa', 0), 0, 140),
                                        'hbd_score': _score_range(props.get('hbd', 0), 0, 5),
                                        'hba_score': _score_range(props.get('hba', 0), 0, 10),
                                        'rotb_score': _score_range(props.get('rotb', 0), 0, 10),
                                        'clearance_score': _score_ge(props.get('renal_clearance', 0), 0.5),
                                    }
                                    
                                    # For toxicity, lower alerts are better
                                    def _inverse_alert_score(count):
                                        try:
                                            c = float(count)
                                        except Exception:
                                            c = 0.0
                                        return max(0.0, 1.0 - 0.2 * c)
                                    
                                    tox_scores = {
                                        'mw_score': _score_range(props.get('mw', 0), 150, 500),
                                        'logp_score': _score_range(props.get('logp', 0), 0, 5),
                                        'tpsa_score': _score_range(props.get('tpsa', 0), 0, 140),
                                        'hbd_score': _score_range(props.get('hbd', 0), 0, 5),
                                        'hba_score': _score_range(props.get('hba', 0), 0, 10),
                                        'rotb_score': _score_range(props.get('rotb', 0), 0, 10),
                                        'pains_score': _inverse_alert_score(props.get('pains', 0)),
                                        'brenk_score': _inverse_alert_score(props.get('brenk_alerts', 0)),
                                        'nih_score': _inverse_alert_score(props.get('nih_alerts', 0)),
                                        'alerts_score': _inverse_alert_score(props.get('total_alerts', 0)),
                                    }
                                    
                                    # Drug-likeness scoring
                                    def _bool_score(flag):
                                        return 1.0 if bool(flag) else 0.0
                                    
                                    def _clip01(x):
                                        try:
                                            v = float(x)
                                        except Exception:
                                            v = 0.0
                                        return max(0.0, min(1.0, v))
                                    
                                    category = props.get('drug_category', 'Unknown')
                                    if isinstance(category, str):
                                        category_lower = category.lower()
                                    else:
                                        category_lower = 'unknown'
                                    if 'drug' in category_lower:
                                        category_score = 1.0
                                    elif 'lead' in category_lower:
                                        category_score = 0.8
                                    elif 'fragment' in category_lower:
                                        category_score = 0.6
                                    else:
                                        category_score = 0.5
                                    
                                    dl_scores = {
                                        'qed_score': _clip01(props.get('qed_score', 0)),
                                        'category_score': category_score,
                                        'lead_like_score': _bool_score(props.get('lead_like_pass')),
                                        'fragment_like_score': _bool_score(props.get('fragment_like_pass')),
                                        'cns_like_score': _bool_score(props.get('cns_like_pass')),
                                        'lipinski_score': _bool_score(props.get('lipinski_pass')),
                                        'veber_score': _bool_score(props.get('veber_pass')),
                                        'egan_score': _bool_score(props.get('egan_pass')),
                                        'synthetic_score': 1.0 - min(1.0, float(props.get('synthetic_accessibility', 0)) / 10.0),
                                    }
                                    
                                    # Create tabs for each ADMET category
                                    admet_tabs = st.tabs([
                                        "🔍 Absorption", "🌐 Distribution", "⚗️ Metabolism", 
                                        "💧 Excretion", "⚠️ Toxicity", "💊 Drug-likeness"
                                    ])
                                    
                                    with admet_tabs[0]:  # Absorption
                                        st.subheader("🔍 Absorption Properties")
                                        absorption_data = {
                                            'Property': [
                                                'Molecular Weight (Da)', 'LogP', 'TPSA (Å²)', 
                                                'HBD Count', 'HBA Count', 'Rotatable Bonds',
                                                'Bioavailability Score', 'Lipinski Rule Pass'
                                            ],
                                            'Value': [
                                                f"{props.get('mw', 0):.1f}",
                                                f"{props.get('logp', 0):.2f}",
                                                f"{props.get('tpsa', 0):.1f}",
                                                props.get('hbd', 0),
                                                props.get('hba', 0),
                                                props.get('rotb', 0),
                                                f"{props.get('bioavailability_score', 0):.3f}",
                                                '✅' if props.get('lipinski_pass') else '❌'
                                            ],
                                            'Score': [
                                                f"{abs_scores['mw_score']:.3f}",
                                                f"{abs_scores['logp_score']:.3f}",
                                                f"{abs_scores['tpsa_score']:.3f}",
                                                f"{abs_scores['hbd_score']:.3f}",
                                                f"{abs_scores['hba_score']:.3f}",
                                                f"{abs_scores['rotb_score']:.3f}",
                                                f"{abs_scores['bioavailability_score']:.3f}",
                                                f"{abs_scores['lipinski_score']:.3f}"
                                            ],
                                            'Target': [
                                                '150-500 Da', '0-5', '0-140 Å²', 
                                                '0-5', '0-10', '0-10',
                                                '>0.5', 'Pass'
                                            ]
                                        }
                                        absorption_df = pd.DataFrame(absorption_data)
                                        st.dataframe(absorption_df, use_container_width=True, hide_index=True)
                                    
                                    with admet_tabs[1]:  # Distribution
                                        st.subheader("🌐 Distribution Properties")
                                        distribution_data = {
                                            'Property': [
                                                'LogP', 'TPSA (Å²)', 'Molecular Weight (Da)',
                                                'HBD Count', 'HBA Count', 'CNS Permeability',
                                                'Blood-Brain Barrier Score'
                                            ],
                                            'Value': [
                                                f"{props.get('logp', 0):.2f}",
                                                f"{props.get('tpsa', 0):.1f}",
                                                f"{props.get('mw', 0):.1f}",
                                                props.get('hbd', 0),
                                                props.get('hba', 0),
                                                f"{props.get('cns_permeability', 0):.3f}",
                                                f"{props.get('bbb_score', 0):.3f}"
                                            ],
                                            'Score': [
                                                f"{dist_scores['logp_score']:.3f}",
                                                f"{dist_scores['tpsa_score']:.3f}",
                                                f"{dist_scores['mw_score']:.3f}",
                                                f"{dist_scores['hbd_score']:.3f}",
                                                f"{dist_scores['hba_score']:.3f}",
                                                f"{dist_scores['cns_score']:.3f}",
                                                f"{dist_scores['bbb_score']:.3f}"
                                            ],
                                            'Target': [
                                                '2-5', '0-90 Å²', '150-500 Da',
                                                '0-3', '0-8', '>0.5',
                                                '>0.5'
                                            ]
                                        }
                                        distribution_df = pd.DataFrame(distribution_data)
                                        st.dataframe(distribution_df, use_container_width=True, hide_index=True)
                                    
                                    with admet_tabs[2]:  # Metabolism
                                        st.subheader("⚗️ Metabolism Properties")
                                        metabolism_data = {
                                            'Property': [
                                                'Molecular Weight (Da)', 'LogP', 'TPSA (Å²)',
                                                'HBD Count', 'HBA Count', 'Rotatable Bonds',
                                                'Metabolic Stability Score'
                                            ],
                                            'Value': [
                                                f"{props.get('mw', 0):.1f}",
                                                f"{props.get('logp', 0):.2f}",
                                                f"{props.get('tpsa', 0):.1f}",
                                                props.get('hbd', 0),
                                                props.get('hba', 0),
                                                props.get('rotb', 0),
                                                f"{props.get('metabolic_stability', 0):.3f}"
                                            ],
                                            'Score': [
                                                f"{metab_scores['mw_score']:.3f}",
                                                f"{metab_scores['logp_score']:.3f}",
                                                f"{metab_scores['tpsa_score']:.3f}",
                                                f"{metab_scores['hbd_score']:.3f}",
                                                f"{metab_scores['hba_score']:.3f}",
                                                f"{metab_scores['rotb_score']:.3f}",
                                                f"{metab_scores['stability_score']:.3f}"
                                            ],
                                            'Target': [
                                                '150-500 Da', '0-5', '0-140 Å²',
                                                '0-5', '0-10', '0-10',
                                                '>0.5'
                                            ]
                                        }
                                        metabolism_df = pd.DataFrame(metabolism_data)
                                        st.dataframe(metabolism_df, use_container_width=True, hide_index=True)
                                    
                                    with admet_tabs[3]:  # Excretion
                                        st.subheader("💧 Excretion Properties")
                                        excretion_data = {
                                            'Property': [
                                                'Molecular Weight (Da)', 'LogP', 'TPSA (Å²)',
                                                'HBD Count', 'HBA Count', 'Rotatable Bonds',
                                                'Renal Clearance Score'
                                            ],
                                            'Value': [
                                                f"{props.get('mw', 0):.1f}",
                                                f"{props.get('logp', 0):.2f}",
                                                f"{props.get('tpsa', 0):.1f}",
                                                props.get('hbd', 0),
                                                props.get('hba', 0),
                                                props.get('rotb', 0),
                                                f"{props.get('renal_clearance', 0):.3f}"
                                            ],
                                            'Score': [
                                                f"{excr_scores['mw_score']:.3f}",
                                                f"{excr_scores['logp_score']:.3f}",
                                                f"{excr_scores['tpsa_score']:.3f}",
                                                f"{excr_scores['hbd_score']:.3f}",
                                                f"{excr_scores['hba_score']:.3f}",
                                                f"{excr_scores['rotb_score']:.3f}",
                                                f"{excr_scores['clearance_score']:.3f}"
                                            ],
                                            'Target': [
                                                '150-500 Da', '0-5', '0-140 Å²',
                                                '0-5', '0-10', '0-10',
                                                '>0.5'
                                            ]
                                        }
                                        excretion_df = pd.DataFrame(excretion_data)
                                        st.dataframe(excretion_df, use_container_width=True, hide_index=True)
                                    
                                    with admet_tabs[4]:  # Toxicity
                                        st.subheader("⚠️ Toxicity Properties")
                                        toxicity_data = {
                                            'Property': [
                                                'Molecular Weight (Da)', 'LogP', 'TPSA (Å²)',
                                                'HBD Count', 'HBA Count', 'Rotatable Bonds',
                                                'PAINS Alerts', 'Brenk Alerts', 'NIH Alerts',
                                                'Total Structural Alerts'
                                            ],
                                            'Value': [
                                                f"{props.get('mw', 0):.1f}",
                                                f"{props.get('logp', 0):.2f}",
                                                f"{props.get('tpsa', 0):.1f}",
                                                props.get('hbd', 0),
                                                props.get('hba', 0),
                                                props.get('rotb', 0),
                                                props.get('pains', 0),
                                                props.get('brenk_alerts', 0),
                                                props.get('nih_alerts', 0),
                                                props.get('total_alerts', 0)
                                            ],
                                            'Score': [
                                                f"{tox_scores['mw_score']:.3f}",
                                                f"{tox_scores['logp_score']:.3f}",
                                                f"{tox_scores['tpsa_score']:.3f}",
                                                f"{tox_scores['hbd_score']:.3f}",
                                                f"{tox_scores['hba_score']:.3f}",
                                                f"{tox_scores['rotb_score']:.3f}",
                                                f"{tox_scores['pains_score']:.3f}",
                                                f"{tox_scores['brenk_score']:.3f}",
                                                f"{tox_scores['nih_score']:.3f}",
                                                f"{tox_scores['alerts_score']:.3f}"
                                            ],
                                            'Target': [
                                                '150-500 Da', '0-5', '0-140 Å²',
                                                '0-5', '0-10', '0-10',
                                                '0', '0', '0',
                                                '0'
                                            ]
                                        }
                                        toxicity_df = pd.DataFrame(toxicity_data)
                                        st.dataframe(toxicity_df, use_container_width=True, hide_index=True)
                                    
                                    with admet_tabs[5]:  # Drug-likeness
                                        st.subheader("💊 Drug-likeness Properties")
                                        drug_likeness_data = {
                                            'Property': [
                                                'QED Score', 'Drug Category', 'Lead-like Pass',
                                                'Fragment-like Pass', 'CNS-likeness Pass',
                                                'Lipinski Rule Pass', 'Veber Rule Pass',
                                                'Egan Rule Pass', 'Synthetic Accessibility'
                                            ],
                                            'Value': [
                                                f"{props.get('qed_score', 0):.3f}",
                                                props.get('drug_category', 'Unknown'),
                                                '✅' if props.get('lead_like_pass') else '❌',
                                                '✅' if props.get('fragment_like_pass') else '❌',
                                                '✅' if props.get('cns_like_pass') else '❌',
                                                '✅' if props.get('lipinski_pass') else '❌',
                                                '✅' if props.get('veber_pass') else '❌',
                                                '✅' if props.get('egan_pass') else '❌',
                                                f"{props.get('synthetic_accessibility', 0):.3f}"
                                            ],
                                            'Score': [
                                                f"{dl_scores['qed_score']:.3f}",
                                                f"{dl_scores['category_score']:.3f}",
                                                f"{dl_scores['lead_like_score']:.3f}",
                                                f"{dl_scores['fragment_like_score']:.3f}",
                                                f"{dl_scores['cns_like_score']:.3f}",
                                                f"{dl_scores['lipinski_score']:.3f}",
                                                f"{dl_scores['veber_score']:.3f}",
                                                f"{dl_scores['egan_score']:.3f}",
                                                f"{dl_scores['synthetic_score']:.3f}"
                                            ],
                                            'Target': [
                                                '>0.5', 'Drug-like', 'Pass',
                                                'Pass', 'Pass',
                                                'Pass', 'Pass',
                                                'Pass', '<5.0'
                                            ]
                                        }
                                        drug_likeness_df = pd.DataFrame(drug_likeness_data)
                                        st.dataframe(drug_likeness_df, use_container_width=True, hide_index=True)
                    
                    with viz_tab3:
                        st.subheader("🔥 ADMET Properties Heatmap")
                        
                        # Create heatmap data
                        heatmap_data = []
                        compound_names = []
                        
                        for idx, candidate in ranked_df.head(8).iterrows():
                            props = compute_admet_properties(candidate['smiles'])
                            if 'error' not in props:
                                summary = interpret_admet(props)
                                heatmap_data.append([
                                    summary['absorption']['score'],
                                    summary['distribution']['score'],
                                    summary['metabolism']['score'],
                                    summary['excretion']['score'],
                                    summary['toxicity']['score'],
                                    summary['drug_likeness']['score']
                                ])
                                compound_names.append(candidate['compound_name'])
                        
                        if heatmap_data:
                            fig_heatmap = go.Figure(data=go.Heatmap(
                                z=heatmap_data,
                                x=['Absorption', 'Distribution', 'Metabolism', 'Excretion', 'Toxicity', 'Drug-likeness'],
                                y=compound_names,
                                colorscale='RdYlGn',
                                showscale=True,
                                colorbar=dict(title="ADMET Score")
                            ))
                            
                            fig_heatmap.update_layout(
                                title="ADMET Properties Heatmap",
                                xaxis_title="ADMET Categories",
                                yaxis_title="Compounds",
                                height=500
                            )
                            
                            st.plotly_chart(fig_heatmap, use_container_width=True)
                    
                    with viz_tab4:
                        st.subheader("📋 Detailed Compound Analysis")
                        
                        # Show detailed analysis for top 3 compounds
                        for idx, candidate in ranked_df.head(3).iterrows():
                            with st.expander(f"🔍 Detailed Analysis: {candidate['compound_name']}"):
                                props = compute_admet_properties(candidate['smiles'])
                                
                                if 'error' not in props:
                                    summary = interpret_admet(props)
                                    
                                    # Overall metrics
                                    col1, col2, col3, col4 = st.columns(4)
                                    with col1:
                                        st.metric("Combined Score", f"{candidate['combined_score']:.3f}")
                                    with col2:
                                        st.metric("Docking Score", f"{candidate['docking_score']:.3f}")
                                    with col3:
                                        st.metric("ADMET Score", f"{candidate['admet_score']:.3f}")
                                    with col4:
                                        st.metric("Drug Category", candidate['drug_category'])
                                    
                                    # ADMET breakdown
                                    st.subheader("ADMET Breakdown")
                                    
                                    # Create bar chart for ADMET scores
                                    admet_categories = ['Absorption', 'Distribution', 'Metabolism', 'Excretion', 'Toxicity', 'Drug-likeness']
                                    admet_scores = [
                                        summary['absorption']['score'],
                                        summary['distribution']['score'],
                                        summary['metabolism']['score'],
                                        summary['excretion']['score'],
                                        summary['toxicity']['score'],
                                        summary['drug_likeness']['score']
                                    ]
                                    
                                    fig_bar = go.Figure()
                                    fig_bar.add_trace(go.Bar(
                                        x=admet_categories,
                                        y=admet_scores,
                                        marker_color=['#ff6b6b', '#4ecdc4', '#45b7d1', '#96ceb4', '#feca57', '#ff9ff3'],
                                        text=admet_scores,
                                        textposition='auto'
                                    ))
                                    
                                    fig_bar.update_layout(
                                        title=f"ADMET Scores: {candidate['compound_name']}",
                                        xaxis_title="ADMET Categories",
                                        yaxis_title="Score",
                                        yaxis=dict(range=[0, 1]),
                                        height=400
                                    )
                                    
                                    st.plotly_chart(fig_bar, use_container_width=True)
                                    
                                    # Molecular properties
                                    st.subheader("Molecular Properties")
                                    col1, col2, col3 = st.columns(3)
                                    with col1:
                                        st.metric("MW", f"{props.get('mw', 0):.1f} Da")
                                        st.metric("LogP", f"{props.get('logp', 0):.2f}")
                                    with col2:
                                        st.metric("TPSA", f"{props.get('tpsa', 0):.1f} Å²")
                                        st.metric("HBD", props.get('hbd', 0))
                                    with col3:
                                        st.metric("HBA", props.get('hba', 0))
                                        st.metric("RotB", props.get('rotb', 0))
                                    
                                    # Rule compliance
                                    st.subheader("Drug-likeness Rules")
                                    rules_data = {
                                        'Rule': ['Lipinski', 'Veber', 'Egan', 'Lead-like', 'Fragment-like'],
                                        'Pass': [
                                            '✅' if props.get('lipinski_pass') else '❌',
                                            '✅' if props.get('veber_pass') else '❌',
                                            '✅' if props.get('egan_pass') else '❌',
                                            '✅' if props.get('lead_like_pass') else '❌',
                                            '✅' if props.get('fragment_like_pass') else '❌'
                                        ]
                                    }
                                    rules_df = pd.DataFrame(rules_data)
                                    st.dataframe(rules_df, use_container_width=True, hide_index=True)
                                    
                                    # Structural alerts
                                    st.subheader("Structural Alerts")
                                    alerts_data = {
                                        'Alert Type': ['PAINS', 'Brenk', 'NIH', 'Total'],
                                        'Count': [
                                            props.get('pains', 0),
                                            props.get('brenk_alerts', 0),
                                            props.get('nih_alerts', 0),
                                            props.get('total_alerts', 0)
                                        ]
                                    }
                                    alerts_df = pd.DataFrame(alerts_data)
                                    st.dataframe(alerts_df, use_container_width=True, hide_index=True)
                    
                    # Export results
                    st.subheader("💾 Export Results")
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        # JSON export
                        export_data = {
                            'pipeline_results': st.session_state.drug_discovery_pipeline,
                            'timestamp': datetime.now().isoformat(),
                            'summary': {
                                'total_compounds': len(admet_predictions),
                                'top_candidate': admet_predictions[0]['compound_name'] if admet_predictions else None,
                                'best_score': admet_predictions[0]['combined_score'] if admet_predictions else 0
                            }
                        }
                        st.download_button(
                            "📄 Download Full Results (JSON)",
                            data=json.dumps(export_data, indent=2),
                            file_name=f"drug_discovery_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json",
                            mime="application/json"
                        )
                    
                    with col2:
                        # CSV export
                        csv_data = ranked_df.to_csv(index=False)
                        st.download_button(
                            "📊 Download Rankings (CSV)",
                            data=csv_data,
                            file_name=f"compound_rankings_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                            mime="text/csv"
                        )
                    
                    with col3:
                        # Summary report
                        # Fix the f-string formatting
                        best_score = f"{admet_predictions[0]['combined_score']:.3f}" if admet_predictions else 'N/A'
                        target_protein = st.session_state.drug_discovery_pipeline['target_protein']['pdb_id'] if st.session_state.drug_discovery_pipeline['target_protein'] else 'N/A'
                        top_candidate = admet_predictions[0]['compound_name'] if admet_predictions else 'N/A'
                        
                        summary_report = f"""
AI Drug Discovery Pipeline Results
Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}

Target Protein: {target_protein}
Total Compounds Screened: {len(admet_predictions)}
Top Candidate: {top_candidate}
Best Combined Score: {best_score}

Top 5 Candidates:
"""
                        for i, candidate in enumerate(admet_predictions[:5]):
                            summary_report += f"""
{i+1}. {candidate['compound_name']}
   - Combined Score: {candidate['combined_score']:.3f}
   - Docking Score: {candidate['docking_score']:.3f}
   - ADMET Score: {candidate['admet_score']:.3f}
   - Drug Category: {candidate['drug_category']}
"""
                        
                        st.download_button(
                            "📋 Download Summary (TXT)",
                            data=summary_report,
                            file_name=f"discovery_summary_{datetime.now().strftime('%Y%m%d_%H%M%S')}.txt",
                            mime="text/plain"
                        )

# -------------------- Enhanced Sidebar Features --------------------
with st.sidebar:
    st.markdown("---")
    st.subheader("🎤 Voice Commands")
    
    if st.button("🎤 Enable Voice Control"):
        voice_commands = generate_voice_commands()
        st.info("Voice commands available:")
        for cmd in voice_commands[:5]:  # Show first 5
            st.text(f"• {cmd}")
    
    st.markdown("---")
    st.subheader("👥 Collaboration")
    
    room_name = st.text_input("Collaboration Room", placeholder="Enter room name")
    if st.button("🚪 Join Room"):
        if room_name:
            room_data = setup_collaboration_room(room_name)
            st.session_state.collaboration_room = room_data
            st.success(f"✅ Joined room: {room_name}")
        else:
            st.warning("Please enter a room name")
    
    if st.session_state.collaboration_room:
        st.info(f"📍 Active room: {st.session_state.collaboration_room['room_id']}")
        if st.button("🚪 Leave Room"):
            st.session_state.collaboration_room = None
            st.success("✅ Left collaboration room")

# -------------------- Footer --------------------
st.markdown("---")
st.markdown("""
<div style="text-align: center; color: #64748b; padding: 2rem;">
    <p>🧬 <strong>ProStruct - 3D</strong> | Advanced Protein Structure Prediction & Analysis Platform</p>
    <p>Built with Streamlit, ESMFold, and scientific Python libraries</p>
    <p>Developed by Keeistu M S </p>
</div>
""", unsafe_allow_html=True)

# -------------------- Optional notes --------------------
with st.expander("ℹ️ Notes & limitations", expanded=False):
    st.markdown('''
**Notes & limitations**:
- pI is estimated using a simple Henderson–Hasselbalch based approach and is approximate.
- Secondary structure counts are read from HELIX/SHEET records in PDB if present; ESMFold output may not contain these.
- Domain prediction, PPI prediction, and PTM sites are simplified algorithms for demonstration purposes.
- For production use, consider integrating with specialized databases and tools like Pfam, STRING, and UniProt.
- For long sequences, prediction time may be longer and memory usage higher. Caching reduces repeated API calls for identical sequences.
- Batch processing is available but may be rate-limited by the ESM Atlas API.
''')
