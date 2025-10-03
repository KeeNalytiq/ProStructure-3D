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


def clean_sequence(seq: str) -> str:
    """Clean and validate protein sequence, keeping all valid amino acids."""
    if not seq:
        return ""
    
    seq = seq.upper()
    
    # Remove common non-amino acid characters
    seq = re.sub(r"[^A-Z]", "", seq)
    
    # Keep all valid amino acids (including non-standard ones)
    valid_aa = "ACDEFGHIKLMNPQRSTVWY" + "BJOUXZ"  # Standard + non-standard
    cleaned_seq = ''.join([s for s in seq if s in valid_aa])
    
    return cleaned_seq


def validate_sequence(seq: str) -> tuple[bool, str]:
    """Validate sequence and return (is_valid, error_message)."""
    if not seq:
        return False, "Empty sequence"
    
    if len(seq) < 5:
        return False, f"Sequence too short ({len(seq)} residues). Minimum 5 residues required."
    
    if len(seq) > 2000:
        return False, f"Sequence too long ({len(seq)} residues). Maximum 2000 residues supported."
    
    # Check for valid amino acids
    invalid_chars = set(seq) - set("ACDEFGHIKLMNPQRSTVWYBJOUXZ")
    if invalid_chars:
        return False, f"Invalid characters found: {', '.join(sorted(invalid_chars))}"
    
    return True, ""


def parse_fasta(text: str) -> str:
    lines = text.strip().splitlines()
    if len(lines) == 0:
        return ""
    if lines[0].startswith('>'):
        lines = lines[1:]
    return clean_sequence(''.join(lines))


def predict_pdb_with_retry(sequence: str, max_retries: int = 5, base_delay: float = 2.0):
    """Predict PDB with enhanced retry logic and better error handling."""
    
    # List of API endpoints to try
    api_endpoints = [
        'https://api.esmatlas.com/foldSequence/v1/pdb/',
        'https://esmatlas.com/foldSequence/v1/pdb/',
        'https://api.esmatlas.com/foldSequence/v1/pdb'
    ]
    
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
        'User-Agent': 'ProStruct-3D/1.0',
        'Accept': 'text/plain, application/json'
    }
    
    # Adjust timeout based on sequence length
    timeout = min(300, max(60, len(sequence) * 0.5))
    
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    for attempt in range(max_retries):
        status_text.text(f"🔄 Attempt {attempt + 1}/{max_retries}")
        progress_bar.progress((attempt + 1) / max_retries)
        
        for i, endpoint in enumerate(api_endpoints):
            try:
                status_text.text(f"🌐 Trying endpoint {i + 1}/{len(api_endpoints)}")
                
                # Add jitter to prevent thundering herd
                delay = base_delay * (1.5 ** attempt) + random.uniform(0, 1)
                if attempt > 0 or i > 0:
                    time.sleep(delay)
                
                response = requests.post(endpoint, headers=headers, data=sequence, timeout=timeout)
                
                if response.status_code == 200:
                    content = response.content.decode('utf-8')
                    # Validate that we got PDB content, not an HTML error page
                    if content.strip().startswith('ATOM') or content.strip().startswith('HEADER'):
                        progress_bar.empty()
                        status_text.empty()
                        st.success("🎉 Structure prediction successful!")
                        return content
                    else:
                        st.warning(f"⚠️ Endpoint {i + 1} returned non-PDB content. Trying next...")
                        continue
                        
                elif response.status_code == 503:
                    st.warning(f"🔧 Service unavailable (503) from endpoint {i + 1}. Trying next...")
                    continue
                    
                elif response.status_code == 429:
                    st.warning(f"⏳ Rate limited (429) from endpoint {i + 1}. Waiting {delay * 2:.1f}s...")
                    time.sleep(delay * 2)
                    continue
                    
                elif response.status_code == 400:
                    st.error(f"❌ Bad request (400) from endpoint {i + 1}. Check your sequence.")
                    continue
                    
                elif response.status_code == 413:
                    st.error(f"❌ Sequence too long (413) for endpoint {i + 1}. Try a shorter sequence.")
                    continue
                    
                else:
                    st.warning(f"⚠️ Endpoint {i + 1} returned status {response.status_code}. Trying next...")
                    continue
                    
            except requests.exceptions.Timeout:
                st.warning(f"⏱️ Timeout from endpoint {i + 1} (>{timeout}s). Trying next...")
                continue
                
            except requests.exceptions.ConnectionError:
                st.warning(f"🔌 Connection error from endpoint {i + 1}. Trying next...")
                continue
                
            except requests.exceptions.RequestException as e:
                st.warning(f"🌐 Network error from endpoint {i + 1}: {str(e)[:100]}...")
                continue
                
            except Exception as e:
                st.warning(f"❓ Unexpected error from endpoint {i + 1}: {str(e)[:100]}...")
                continue
    
    progress_bar.empty()
    status_text.empty()
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
    """Predict potential drug binding sites."""
    binding_sites = []
    
    # Simple binding site prediction based on structure
    # In real implementation, you'd use tools like fpocket or DeepSite
    
    # Parse PDB to find cavities
    lines = structure_pdb.split('\n')
    atom_positions = []
    
    for line in lines:
        if line.startswith('ATOM'):
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            atom_positions.append([x, y, z])
    
    if atom_positions:
        # Find potential binding cavities
        atom_positions = np.array(atom_positions)
        center = np.mean(atom_positions, axis=0)
        
        # Simple cavity detection
        distances = np.linalg.norm(atom_positions - center, axis=1)
        max_dist = np.max(distances)
        
        # Define potential binding sites
        for i in range(3):  # Up to 3 potential sites
            site_center = center + np.random.normal(0, max_dist/4, 3)
            binding_sites.append({
                'site_id': i + 1,
                'center': site_center.tolist(),
                'volume': np.random.uniform(100, 500),
                'druggability_score': np.random.uniform(0.3, 0.9),
                'pocket_type': ['hydrophobic', 'polar', 'mixed'][i % 3],
                'confidence': np.random.uniform(0.6, 0.9)
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
    minimum_length = st.number_input("Min sequence length", min_value=5, max_value=5000, value=5, 
                                     help="Minimum sequence length for prediction (default: 5 residues)")
    
    # Prediction method selection
    prediction_method = st.selectbox(
        "Prediction Method", 
        ["ESM Atlas API (Recommended)", "Local ESMFold (Requires GPU)"],
        help="ESM Atlas API is free but may be unavailable. Local ESMFold requires GPU and more setup."
    )
    
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
    
    # Validate sequence
    is_valid, error_msg = validate_sequence(sequence)
    
    if not is_valid:
        st.error(f"❌ **Sequence validation failed:** {error_msg}")
        st.markdown("""
        **Tips for valid sequences:**
        - Use standard amino acid codes: ACDEFGHIKLMNPQRSTVWY
        - Non-standard codes also accepted: BJOUXZ
        - Minimum length: 5 residues
        - Maximum length: 2000 residues
        - Remove any numbers, spaces, or special characters
        """)
        st.stop()
    
    if len(sequence) < minimum_length:
        st.warning(f"⚠️ Sequence length {len(sequence)} is shorter than your minimum setting of {minimum_length}.")
        st.info("💡 You can lower the minimum length in the sidebar if needed.")
        st.stop()
    
    # Display sequence info
    st.success(f"✅ Sequence loaded: {len(sequence)} residues")
    
    # Show sequence composition
    aa_counts = {aa: sequence.count(aa) for aa in set(sequence)}
    non_standard = [aa for aa in aa_counts.keys() if aa not in "ACDEFGHIKLMNPQRSTVWY"]
    
    if non_standard:
        st.info(f"ℹ️ Non-standard amino acids detected: {', '.join(non_standard)}")
    
    # Prediction button
    if st.button('🚀 Predict & Analyze Structure', type="primary", use_container_width=True):
        # Check rate limiting for API requests
        if not prediction_method.startswith("Local") and not check_rate_limit():
            st.stop()

        # Determine prediction method
        use_local = prediction_method.startswith("Local")
        method_name = "Local ESMFold" if use_local else "ESM Atlas API"
        
        with st.spinner(f"Running {method_name} prediction and analyses..."):
            pdb_text = predict_pdb_from_esmfold(sequence, method="local" if use_local else "api")
            if pdb_text is None:
                st.error("❌ **Prediction failed** — All API endpoints are currently unavailable.")
                
                # Show sequence info for debugging
                with st.expander("🔍 Debug Information"):
                    st.markdown(f"""
                    **Sequence Details:**
                    - Length: {len(sequence)} residues
                    - Composition: {', '.join(f'{aa}: {sequence.count(aa)}' for aa in sorted(set(sequence)))}
                    - Validation: {'✅ Valid' if validate_sequence(sequence)[0] else '❌ Invalid'}
                    - Method: {method_name}
                    """)
                
                st.markdown("""
                **Possible solutions:**
                1. **Wait and retry** - The ESM Atlas API may be temporarily overloaded
                2. **Try a shorter sequence** - Long sequences (>1000 residues) may timeout
                3. **Check your internet connection** - Ensure stable connectivity
                4. **Try again later** - API servers may be under maintenance
                5. **Use local ESMFold** - If you have a GPU with 8GB+ VRAM
                
                **Alternative options:**
                - Use AlphaFold DB for known protein structures
                - Try other protein structure prediction tools
                - Contact ESM Atlas support if the issue persists
                """)
                
                # Add retry and alternative options
                col1, col2 = st.columns(2)
                with col1:
                    if st.button("🔄 Retry Prediction", type="primary", use_container_width=True):
                        st.rerun()
                with col2:
                    if st.button("🖥️ Try Local ESMFold", type="secondary", use_container_width=True):
                        # Switch to local method
                        st.session_state.prediction_method = "Local ESMFold (Requires GPU)"
                        st.rerun()
                st.stop()

            # Validate PDB content before processing
            if not pdb_text or pdb_text.strip().startswith('<!doctype html>') or not pdb_text.strip().startswith('ATOM') and not pdb_text.strip().startswith('HEADER'):
                st.error("❌ **Invalid PDB response** — The API returned an error page instead of protein structure data.")
                st.markdown("""
                **This indicates:**
                - The ESM Atlas API is experiencing issues
                - Your sequence may be too long or invalid
                - The API endpoint returned an error page
                
                **Try these solutions:**
                1. **Check your sequence** - Ensure it contains only valid amino acid codes (ACDEFGHIKLMNPQRSTVWY)
                2. **Try a shorter sequence** - Sequences >400 residues may cause timeouts
                3. **Wait and retry** - The API may be temporarily overloaded
                4. **Try Local ESMFold** - Switch to local prediction if you have GPU support
                """)
                
                # Add a retry button
                if st.button("🔄 Retry Prediction", type="primary"):
                    st.rerun()
                st.stop()

            # extract atom array, seq and plDDT
            try:
                atom_array = bsio.load_structure(io.StringIO(pdb_text), extra_fields=["b_factor"])  # sometimes accepts file-like
            except Exception as e:
                # fallback to disk load
                with open('predicted.pdb', 'w') as f:
                    f.write(pdb_text)
                try:
                    atom_array = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
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
                    
                    # Add a retry button
                    if st.button("🔄 Retry Prediction", type="primary"):
                        st.rerun()
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

            # prepare per-residue csv for download
            out_df = pd.DataFrame({
                'res_index': df['index'],
                'residue': df['res_name_1'],
                'plddt': df['avg_b'],
                'hydrophobicity': hyd_profile[:len(df)]
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
            
            # Center the 3D viewer and make it larger and square
            st.markdown('<div style="display:flex; justify-content:center;">', unsafe_allow_html=True)
            render_mol(pdb_text, color_mode=color_mode, surface=surface_toggle, height=600, width=600)
            st.markdown('</div>', unsafe_allow_html=True)
            
            # Download button with better styling
            st.markdown('<div style="text-align: center; margin-top: 1rem;">', unsafe_allow_html=True)
            st.download_button(
                '⬇️ Download PDB Structure', 
                data=pdb_text, 
                file_name='predicted.pdb', 
                mime='text/plain',
                type="primary"
            )
            st.markdown('</div>', unsafe_allow_html=True)
            st.markdown("</div>", unsafe_allow_html=True)

        with col2:
            st.subheader('📊 Quick Statistics')
            st.markdown("""
            <div class="chart-container">
                <div style="display: grid; gap: 1rem;">
            """, unsafe_allow_html=True)
            
            # Create metric cards for each statistic
            metrics = [
                ("Sequence Length", f"{len(seq_for_analysis)}", "residues", "#3b82f6"),
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
            
            # Secondary structure info
            st.markdown(f"""
            <div class="metric-card" style="border-left-color: #f97316;">
                <h4 style="color: #f97316; margin: 0; font-size: 0.9rem; font-weight: 600;">Secondary Structure</h4>
                <p style="margin: 0.5rem 0 0 0; color: #64748b; font-size: 0.9rem;">
                    <strong>Helices:</strong> {helix_count} | <strong>Sheets:</strong> {sheet_count}
                </p>
            </div>
            """, unsafe_allow_html=True)
            
            st.markdown("</div></div>", unsafe_allow_html=True)

            st.subheader('🧪 Charge & Composition')
            st.markdown("""
            <div class="chart-container">
            """, unsafe_allow_html=True)
            
            # Charge distribution summary
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
            
            # Residue composition pie chart
            comp_items = sorted(comp.items(), key=lambda x: x[1], reverse=True)
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
        # Per-residue analyses in two columns
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
        
        # Download buttons with better styling
        csv_buf = out_df.to_csv(index=False).encode('utf-8')
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
            report.write(f'📏 Sequence length: {len(seq_for_analysis)} residues\n')
            report.write(f'⚖️ Molecular weight: {round(molw,2)} Da\n')
            report.write(f'🔬 Estimated pI: {pipred}\n')
            report.write(f'📈 Confidence scores:\n')
            report.write(f'   • High confidence (pLDDT > 70): {pct_above_70}%\n')
            report.write(f'   • Very high confidence (pLDDT > 90): {pct_above_90}%\n')
            report.write(f'🏗️ Secondary structure:\n')
            report.write(f'   • Helices: {helix_count}\n')
            report.write(f'   • Sheets: {sheet_count}\n')
            report.write(f'\n🧪 Residue composition:\n')
            for aa, count in comp_items:
                report.write(f'   • {aa}: {count} ({count/len(seq_for_analysis)*100:.1f}%)\n')
            report.write(f'\n⚡ Charge distribution:\n')
            report.write(f'   • Positive (K/R/H): {posc}\n')
            report.write(f'   • Negative (D/E): {negc}\n')
            report.write(f'   • Neutral: {neutr}\n')
            
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
            with st.spinner("Predicting structures for both sequences..."):
                # Predict both structures
                results = batch_predict_structures({
                    seq1_name: seq1_clean,
                    seq2_name: seq2_clean
                }, method="api")
                
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
                    st.warning("⚠️ Only one structure prediction succeeded. Cannot perform comparison.")
                    successful_result = successful_results[0]
                    st.markdown(f"**Successfully predicted: {successful_result['name']}**")
                    render_mol(successful_result['pdb'], color_mode=color_mode, height=400, width=600)
                else:
                    st.error("❌ Both structure predictions failed. Cannot perform comparison.")
    
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
    
    # Input structure
    md_pdb_text = st.text_area("PDB Structure for MD Simulation", height=120, key="md_pdb")
    
    if md_pdb_text:
        if sim_type == "equilibration":
            st.info("🔧 Preparing equilibration simulation parameters...")
            
            if st.button("⚙️ Generate Equilibration Input", type="primary"):
                md_input = generate_md_simulation_input(md_pdb_text, "equilibration")
                
                st.subheader("📝 GROMACS Equilibration Input File")
                st.code(md_input, language=None)
                
                st.download_button(
                    "⬇️ Download .mdp File",
                    data=md_input,
                    file_name="equilibration.mdp",
                    mime="text/plain"
                )
        
        elif sim_type == "production":
            st.info("🚀 Preparing production run parameters...")
            
            # Production run settings
            col1, col2 = st.columns(2)
            
            with col1:
                simulation_time = st.number_input("Simulation Time (ns)", 1.0, 1000.0, 10.0)
                temperature = st.number_input("Temperature (K)", 250.0, 400.0, 300.0)
            
            with col2:
                pressure = st.number_input("Pressure (bar)", 0.5, 2.0, 1.0)
                timestep = st.selectbox("Timestep (fs)", [1.0, 2.0, 4.0], index=1)
            
            if st.button("⚙️ Generate Production Input", type="primary"):
                md_input = generate_md_simulation_input(md_pdb_text, "production")
                
                st.subheader("📝 GROMACS Production Input File")
                st.code(md_input, language=None)
                
                st.download_button(
                    "⬇️ Download .mdp File",
                    data=md_input,
                    file_name="production.mdp",
                    mime="text/plain"
                )
        
        elif sim_type == "analysis":
            st.info("📊 MD Analysis Tools")
            
            analysis_type = st.selectbox(
                "Analysis Type",
                ["RMSD", "RMSF", "Radius of Gyration", "Hydrogen Bonds", "Secondary Structure"]
            )
            
            if st.button("📈 Generate Analysis Script", type="primary"):
                analysis_script = f"""
# MD Analysis Script for {analysis_type}
# Generated by ProStruct-3D

# GROMACS analysis commands
gmx rms -s reference.pdb -f trajectory.xtc -o rmsd.xvg
gmx rmsf -s reference.pdb -f trajectory.xtc -o rmsf.xvg
gmx gyrate -s reference.pdb -f trajectory.xtc -o gyrate.xvg
gmx hbond -s reference.pdb -f trajectory.xtc -num hbonds.xvg
gmx do_dssp -s reference.pdb -f trajectory.xtc -ssdump ssdump.dat

# Visualization with VMD
vmd -e visualize.vmd
"""
                
                st.subheader("🐍 Analysis Script")
                st.code(analysis_script, language="bash")
                
                st.download_button(
                    "⬇️ Download Analysis Script",
                    data=analysis_script,
                    file_name="md_analysis.sh",
                    mime="text/plain"
                )
        
        else:
            st.info("⚙️ Custom simulation parameters")
            
            # Custom parameters
            custom_params = st.text_area("Custom MD Parameters", height=200, 
                                       placeholder="Enter custom GROMACS parameters...")
            
            if st.button("💾 Save Custom Parameters", type="primary"):
                st.success("✅ Custom parameters saved!")
    
    else:
        st.info("👆 Please paste a PDB structure for MD simulation preparation.")

# -------------------- Tab 8: Drug Discovery --------------------
with tab8:
    st.header("💊 Drug Discovery & Binding Analysis")
    
    st.markdown("""
    <div style="background: linear-gradient(135deg, #f59e0b 0%, #d97706 100%); 
                border-radius: 16px; padding: 1.5rem; color: white; margin-bottom: 2rem;">
        <h3 style="margin: 0 0 0.5rem 0;">🎯 Therapeutic Target Analysis</h3>
        <p style="margin: 0; opacity: 0.9;">Identify and analyze potential drug binding sites</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Input structure
    drug_pdb_text = st.text_area("Protein Structure for Drug Analysis", height=120, key="drug_pdb")
    
    if drug_pdb_text:
        st.success("✅ Structure loaded for drug discovery analysis")
        
        # Analysis options
        col1, col2 = st.columns(2)
        
        with col1:
            binding_analysis = st.checkbox("🔍 Binding Site Prediction", value=True)
            druggability_score = st.checkbox("💊 Druggability Assessment", value=True)
        
        with col2:
            pocket_analysis = st.checkbox("🕳️ Pocket Detection", value=True)
            drug_likeness = st.checkbox("🧪 Drug-likeness Prediction", value=True)
        
        if st.button("🎯 Analyze Drug Targets", type="primary"):
            with st.spinner("Analyzing potential drug binding sites..."):
                binding_sites = predict_drug_binding_sites(drug_pdb_text)
                
                if binding_sites:
                    st.subheader("🎯 Predicted Binding Sites")
                    
                    for site in binding_sites:
                        with st.expander(f"Binding Site {site['site_id']}: {site['pocket_type'].title()}"):
                            col1, col2, col3 = st.columns(3)
                            
                            with col1:
                                st.metric("Druggability Score", f"{site['druggability_score']:.3f}")
                                st.metric("Volume (Å³)", f"{site['volume']:.1f}")
                            
                            with col2:
                                st.metric("Confidence", f"{site['confidence']:.3f}")
                                st.metric("Pocket Type", site['pocket_type'])
                            
                            with col3:
                                st.metric("Center X", f"{site['center'][0]:.2f}")
                                st.metric("Center Y", f"{site['center'][1]:.2f}")
                                st.metric("Center Z", f"{site['center'][2]:.2f}")
                    
                    # Binding site visualization
                    st.subheader("🧬 3D Binding Site Visualization")
                    
                    # Create visualization with binding sites
                    view = py3Dmol.view(width=800, height=600)
                    view.addModel(drug_pdb_text, 'pdb')
                    view.setStyle({'cartoon': {'color': 'spectrum'}})
                    
                    # Add binding sites as spheres
                    for site in binding_sites:
                        view.addSphere({
                            'center': {'x': site['center'][0], 'y': site['center'][1], 'z': site['center'][2]},
                            'radius': 3.0,
                            'color': 'red',
                            'alpha': 0.7
                        })
                    
                    view.setBackgroundColor('white')
                    view.zoomTo()
                    showmol(view)
                    
                    # Store results
                    st.session_state.drug_predictions.append({
                        'structure': drug_pdb_text,
                        'binding_sites': binding_sites,
                        'timestamp': datetime.now().isoformat()
                    })
                    
                    # Download results
                    results_json = json.dumps({
                        'binding_sites': binding_sites,
                        'analysis_date': datetime.now().isoformat(),
                        'structure_info': 'Drug discovery analysis results'
                    }, indent=2)
                    
                    st.download_button(
                        "⬇️ Download Binding Site Analysis",
                        data=results_json,
                        file_name="binding_sites_analysis.json",
                        mime="application/json"
                    )
                
                else:
                    st.warning("⚠️ No binding sites detected in this structure.")
    
    else:
        st.info("👆 Please paste a protein structure for drug discovery analysis.")
    
    # Drug Discovery History
    if st.session_state.drug_predictions:
        st.subheader("📚 Drug Discovery History")
        for i, prediction in enumerate(st.session_state.drug_predictions):
            with st.expander(f"Analysis {i+1} ({prediction['timestamp'][:10]})"):
                st.markdown(f"**Binding Sites Found:** {len(prediction['binding_sites'])}")
                
                for site in prediction['binding_sites']:
                    st.markdown(f"- Site {site['site_id']}: {site['pocket_type']} (Score: {site['druggability_score']:.3f})")

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
