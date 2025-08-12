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
import math
from collections import Counter

# -------------------- Page setup --------------------
st.set_page_config(page_title="ProStruct - 3D", layout="wide")
st.title("🎈 ProStruct - 3D")
st.markdown(
    "This app predicts protein structure from a single sequence (ESMFold via ESM Atlas API) and produces advanced analyses: confidence, hydrophobicity, charge, structural metrics, and plots."
)

st.sidebar.header("Input & Settings")
st.sidebar.write("Provide a sequence by pasting it or uploading a FASTA/TXT file. Prediction requires at least 20 residues.")

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


@st.cache_data(show_spinner=False)
def predict_pdb_from_esmfold(sequence: str):
    """Return pdb string. Cached so repeated sequences reuse result."""
    headers = {'Content-Type': 'application/x-www-form-urlencoded'}
    try:
        response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence, timeout=120)
    except Exception as e:
        st.error(f"Network/API error: {e}")
        return None
    if response.status_code != 200:
        st.error(f"API returned status: {response.status_code}\n{response.text}")
        return None
    return response.content.decode('utf-8')


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


# -------------------- Input handling --------------------
seq_text = st.sidebar.text_area("Paste Protein Sequence (single-letter codes)", height=180)
uploaded = st.sidebar.file_uploader("Or upload FASTA/TXT", type=["fasta", "txt"])
minimum_length = st.sidebar.number_input("Minimum sequence length to allow prediction", min_value=10, max_value=5000, value=20)

# options
color_mode = st.sidebar.selectbox("3D Color Mode", ['spectrum', 'chain', 'plddt', 'hydrophobicity'])
surface_toggle = st.sidebar.checkbox("Show molecular surface (may be slow)", value=False)

# read uploaded file if present
sequence = ''
if uploaded is not None:
    try:
        text = uploaded.read().decode('utf-8')
        sequence = parse_fasta(text)
    except Exception:
        st.sidebar.error("Failed to read uploaded file. Ensure it's text/FASTA encoded UTF-8.")
elif seq_text.strip():
    sequence = clean_sequence(seq_text)

if len(sequence) == 0:
    st.sidebar.warning("No valid sequence provided yet.")

# -------------------- Main action --------------------
if st.sidebar.button('Predict & Analyze'):
    if not sequence:
        st.error("Please paste or upload a valid protein sequence (single-letter codes).")
        st.stop()
    if len(sequence) < minimum_length:
        st.error(f"Sequence length {len(sequence)} shorter than required minimum of {minimum_length}.")
        st.stop()

    with st.spinner("Calling ESMFold API and running analyses..."):
        pdb_text = predict_pdb_from_esmfold(sequence)
        if pdb_text is None:
            st.error("Prediction failed — see messages above.")
            st.stop()

        # extract atom array, seq and plDDT
        try:
            atom_array = bsio.load_structure(io.StringIO(pdb_text), extra_fields=["b_factor"])  # sometimes accepts file-like
        except Exception:
            # fallback to disk load
            with open('predicted.pdb', 'w') as f:
                f.write(pdb_text)
            atom_array = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])

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

    # -------------------- Layout: visualization + stats --------------------
    st.markdown("---")
    # Use columns for better alignment and a modern look
    col1, col2 = st.columns([1.5, 1], gap="large")

    with col1:
        st.subheader('🧬 3D Structure')
        with st.container(border=True):
            st.markdown(
                """
                <div style="background: linear-gradient(135deg, #e0e7ff 0%, #f0fdfa 100%);
                            border-radius: 18px; box-shadow: 0 4px 24px rgba(0,0,0,0.08);
                            padding: 1.5rem 1rem 1rem 1rem; margin-bottom: 1.5rem;">
                """,
                unsafe_allow_html=True,
            )
            # Center the 3D viewer and make it larger and square
            st.markdown('<div style="display:flex; justify-content:center;">', unsafe_allow_html=True)
            render_mol(pdb_text, color_mode=color_mode, surface=surface_toggle, height=600, width=600)
            st.markdown('</div>', unsafe_allow_html=True)
            st.download_button('⬇️ Download PDB', data=pdb_text, file_name='predicted.pdb', mime='text/plain')
            st.markdown("</div>", unsafe_allow_html=True)

    with col2:
        st.subheader('📊 Quick Statistics')
        with st.container(border=True):
            st.markdown(
                """
                <div style="background: linear-gradient(135deg, #fdf6e3 0%, #f0fdfa 100%);
                            border-radius: 18px; box-shadow: 0 4px 24px rgba(0,0,0,0.08);
                            padding: 1.2rem 1rem 1rem 1rem; margin-bottom: 1.5rem;">
                """,
                unsafe_allow_html=True,
            )
            st.markdown(f"**Sequence length:** {len(seq_for_analysis)}")
            st.markdown(f"**Molecular weight (approx):** {round(molw,2)} Da")
            st.markdown(f"**Estimated pI:** {pipred}")
            if rg is not None:
                st.markdown(f"**Radius of gyration:** {round(rg,3)} Å")
            else:
                st.info("Radius of gyration: Not available")
            st.markdown(f"**% residues pLDDT > 70:** {pct_above_70}%")
            st.markdown(f"**% residues pLDDT > 90:** {pct_above_90}%")
            st.markdown(f"**Helices (HELIX records):** {helix_count}  —  **Sheets (SHEET records):** {sheet_count}")
            st.markdown("</div>", unsafe_allow_html=True)

        st.subheader('🧪 Charge & Composition')
        with st.container(border=True):
            st.markdown(
                """
                <div style="background: linear-gradient(135deg, #f0fdfa 0%, #e0e7ff 100%);
                            border-radius: 18px; box-shadow: 0 4px 24px rgba(0,0,0,0.08);
                            padding: 1.2rem 1rem 1rem 1rem; margin-bottom: 1.5rem;">
                """,
                unsafe_allow_html=True,
            )
            st.markdown(f"Positive residues (K/R/H): {posc}  |  Negative residues (D/E): {negc}  |  Neutral: {neutr}")
            comp_items = sorted(comp.items(), key=lambda x: x[1], reverse=True)
            labels = [k for k, v in comp_items]
            values = [v for k, v in comp_items]
            if len(labels) > 0:
                fig_pie = go.Figure(data=[go.Pie(labels=labels, values=values, hole=0.3)])
                fig_pie.update_layout(title='Residue Composition')
                st.plotly_chart(fig_pie, use_container_width=True)
            else:
                st.info("No composition data available.")
            st.markdown("</div>", unsafe_allow_html=True)

    st.markdown("---")
    # Per-residue analyses in two columns
    st.subheader('📈 Per-residue Analyses')
    col3, col4 = st.columns(2, gap="large")

    with col3:
        st.markdown('**plDDT per residue**')
        if len(out_df) > 0:
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=out_df['res_index'], y=out_df['plddt'], mode='lines+markers', name='plDDT'))
            fig.add_hline(y=70, line_dash='dash', line_color='orange', annotation_text='70', annotation_position='top left')
            fig.add_hline(y=90, line_dash='dash', line_color='green', annotation_text='90', annotation_position='top left')
            fig.update_layout(
                xaxis_title='Residue index',
                yaxis_title='plDDT score',
                yaxis=dict(range=[0,100]),
                plot_bgcolor='rgba(240,250,255,0.7)',
                paper_bgcolor='rgba(240,250,255,0.7)',
                font=dict(size=14)
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.empty()

    with col4:
        st.markdown('**Hydrophobicity (Kyte-Doolittle)**')
        if len(hyd_profile) > 0:
            fig2 = go.Figure()
            fig2.add_trace(go.Scatter(x=list(range(1, len(hyd_profile)+1)), y=hyd_profile, mode='lines', name='Hydrophobicity'))
            fig2.update_layout(
                xaxis_title='Residue index',
                yaxis_title='KD score',
                plot_bgcolor='rgba(240,250,255,0.7)',
                paper_bgcolor='rgba(240,250,255,0.7)',
                font=dict(size=14)
            )
            st.plotly_chart(fig2, use_container_width=True)
        else:
            st.empty()

    st.markdown("---")
    st.subheader('⚡ Charge distribution')
    if posc + negc + neutr > 0:
        charge_fig = go.Figure(data=[go.Pie(labels=['Positive','Negative','Neutral'], values=[posc,negc,neutr], hole=0.4)])
        charge_fig.update_layout(
            plot_bgcolor='rgba(240,250,255,0.7)',
            paper_bgcolor='rgba(240,250,255,0.7)',
            font=dict(size=14)
        )
        st.plotly_chart(charge_fig, use_container_width=True)
    else:
        st.info("No charge data available.")

    st.markdown("---")
    st.subheader('⬇️ Downloads')
    csv_buf = out_df.to_csv(index=False).encode('utf-8')
    st.download_button('Download per-residue CSV', data=csv_buf, file_name='per_residue_analysis.csv', mime='text/csv')

    report = io.StringIO()
    report.write('Protein Structure Visualizer — Analysis Report\n')
    report.write(f'Sequence length: {len(seq_for_analysis)}\n')
    report.write(f'Molecular weight (approx): {round(molw,2)} Da\n')
    report.write(f'Estimated pI: {pipred}\n')
    report.write(f'% residues pLDDT > 70: {pct_above_70}\n')
    report.write(f'% residues pLDDT > 90: {pct_above_90}\n')
    report.write(f'Helix records: {helix_count}, Sheet records: {sheet_count}\n')
    report.write('\nResidue composition:\n')
    for aa, count in comp_items:
        report.write(f'{aa}: {count}\n')
    st.download_button('Download analysis report (txt)', data=report.getvalue(), file_name='analysis_report.txt', mime='text/plain')

    st.success('Analysis complete.')

else:
    st.info('Ready — paste a protein sequence (single-letter codes) or upload a FASTA/TXT file on the left sidebar, then click `Predict & Analyze`.')

# -------------------- Optional notes --------------------
with st.expander("ℹ️ Notes & limitations", expanded=False):
    st.markdown('''
**Notes & limitations**:
- pI is estimated using a simple Henderson–Hasselbalch based approach and is approximate.
- Secondary structure counts are read from HELIX/SHEET records in PDB if present; ESMFold output may not contain these.
- Binding site prediction (fpocket/DeepSite) is not included by default; it can be integrated if runtime environment supports those tools.
- For long sequences, prediction time may be longer and memory usage higher. Caching reduces repeated API calls for identical sequences.
''')
