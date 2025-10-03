# 🧬 ProStruct - 3D | Advanced Protein Structure Analysis Platform

[![Streamlit](https://img.shields.io/badge/Streamlit-FF4B4B?style=for-the-badge&logo=Streamlit&logoColor=white)](https://share.streamlit.io)
[![Python](https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white)](https://python.org)
[![ESMFold](https://img.shields.io/badge/ESMFold-Meta-FF6B6B?style=for-the-badge)](https://esmatlas.com)
[![License](https://img.shields.io/badge/License-MIT-green.svg?style=for-the-badge)](LICENSE)

> **The Ultimate Protein Structure Prediction & Analysis Platform** - From basic structure prediction to advanced drug discovery and AI-powered protein design.

## 🌟 Overview

**ProStruct - 3D** is a comprehensive, next-generation protein structure analysis platform that combines cutting-edge AI technology with intuitive design. Built on ESMFold's state-of-the-art protein structure prediction, this platform offers everything from basic structure visualization to advanced drug discovery tools.

### 🚀 **Live Demo**
[![Deploy on Streamlit](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://prostruct-3d.streamlit.app/)

## ✨ Key Features

### 🔬 **Core Structure Prediction**
- **ESMFold Integration**: State-of-the-art protein structure prediction
- **Multiple Prediction Methods**: API-based and local GPU predictions
- **Robust Error Handling**: Automatic retry mechanisms and fallback systems
- **Real-time Visualization**: Interactive 3D molecular viewers

### 🤖 **AI-Powered Protein Design**
- **Intelligent Mutation Suggestions**: AI-driven recommendations for stability, activity, and solubility
- **Confidence Scoring**: Advanced algorithms to evaluate mutation impact
- **Design Optimization**: Target-specific protein engineering
- **Design History Tracking**: Complete audit trail of modifications

### 🔄 **Advanced Analysis Tools**
- **Protein Comparison**: Sequence alignment and structural RMSD calculations
- **Domain Analysis**: Pfam domain prediction and visualization
- **PPI Prediction**: Protein-protein interaction networks
- **Evolutionary Conservation**: Conservation analysis and phylogenetic insights

### 🧬 **Molecular Dynamics Simulation**
- **GROMACS Integration**: Complete MD simulation setup
- **Equilibration & Production**: Automated parameter generation
- **Analysis Scripts**: RMSD, RMSF, hydrogen bond analysis
- **Custom Parameters**: Flexible simulation configuration

### 💊 **Drug Discovery Platform**
- **Binding Site Prediction**: Advanced cavity detection algorithms
- **Druggability Assessment**: Comprehensive scoring systems
- **Therapeutic Target Analysis**: Drug-likeness predictions
- **3D Visualization**: Interactive binding site exploration

### 🎨 **Professional UI/UX**
- **Modern Design**: Gradient backgrounds and responsive layout
- **Interactive Charts**: Plotly-powered data visualizations
- **Tabbed Navigation**: Organized workflow management
- **Mobile Responsive**: Optimized for all devices

### 👥 **Collaboration Features**
- **Real-time Sharing**: Collaborative workspaces
- **Voice Commands**: Hands-free operation
- **Export Capabilities**: Multiple format support
- **Session Management**: Persistent data storage

## 🛠️ Installation & Setup

### **Prerequisites**
- Python 3.8+
- 4GB+ RAM (8GB+ recommended)
- GPU with 8GB+ VRAM (for local predictions)

### **Quick Installation**

```bash
# Clone the repository
git clone https://github.com/KeeNalytiq/ProStructure-3D.git
cd ProStructure-3D

# Install dependencies
pip install -r requirements.txt

# Run the application
streamlit run streamlit_app.py
```

### **Local ESMFold Setup (Optional)**

For offline predictions or when API is unavailable:

```bash
# Install PyTorch with CUDA support
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# Install additional dependencies
pip install transformers

# Run with local prediction enabled
streamlit run streamlit_app.py
```

## 🎯 Usage Guide

### **1. Structure Prediction**
- Enter protein sequence or upload FASTA file
- Choose prediction method (API or Local)
- Visualize 3D structure with multiple coloring options
- Analyze confidence scores and structural metrics

### **2. AI Protein Design**
- Input target protein sequence
- Select design objective (stability, activity, solubility)
- Review AI-generated mutation suggestions
- Apply mutations and track design history

### **3. Drug Discovery**
- Upload protein structure (PDB format)
- Run binding site prediction
- Analyze druggability scores
- Export results for further analysis

### **4. Molecular Dynamics**
- Prepare simulation input files
- Configure equilibration and production runs
- Generate analysis scripts
- Download GROMACS-compatible files

## 📊 Technical Specifications

### **Supported Formats**
- **Input**: FASTA, PDB, PDBQT
- **Output**: PDB, JSON, CSV, ZIP
- **Visualization**: PNG, SVG, PyMOL, ChimeraX

### **API Integration**
- **ESM Atlas API**: Free protein structure prediction
- **Rate Limiting**: 1 request per 10 seconds
- **Automatic Retry**: Exponential backoff strategy
- **Multiple Endpoints**: Fallback system for reliability

### **Performance Optimization**
- **Caching**: Intelligent data caching with `@st.cache_data`
- **Lazy Loading**: On-demand library imports
- **Memory Management**: Efficient session state handling
- **Progress Tracking**: Real-time operation feedback

## 🔧 Configuration

### **Environment Variables**
```bash
STREAMLIT_SERVER_PORT=8501
STREAMLIT_SERVER_ADDRESS=0.0.0.0
STREAMLIT_BROWSER_GATHER_USAGE_STATS=false
```

### **Custom Settings**
- Modify `.streamlit/config.toml` for theme customization
- Adjust rate limiting in `streamlit_app.py`
- Configure API endpoints as needed

## 🚀 Deployment Options

### **Streamlit Community Cloud (Recommended)**
1. Push code to GitHub
2. Visit [share.streamlit.io](https://share.streamlit.io)
3. Connect GitHub repository
4. Deploy with one click

### **Other Platforms**
- **Railway**: `railway.json` configuration included
- **Render**: `render.yaml` configuration included
- **Heroku**: `Procfile` configuration included

## 📈 Performance Metrics

- **Prediction Speed**: 30-60 seconds per structure
- **Memory Usage**: 2-4GB RAM typical
- **GPU Requirements**: 8GB+ VRAM for local predictions
- **API Limits**: 1 request per 10 seconds
- **Supported Sequences**: Up to 1024 residues

## 🔬 Scientific Applications

### **Research Use Cases**
- Protein engineering and design
- Drug discovery and development
- Structural biology research
- Bioinformatics education
- Molecular dynamics studies

### **Industry Applications**
- Pharmaceutical research
- Biotechnology development
- Food science applications
- Environmental biotechnology
- Agricultural research

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### **Development Setup**
```bash
# Fork the repository
git clone https://github.com/your-username/ProStructure-3D.git
cd ProStructure-3D

# Install development dependencies
pip install -r requirements-dev.txt

# Run tests
python -m pytest tests/

# Start development server
streamlit run streamlit_app.py
```

## 📚 Documentation

- **API Reference**: [docs/api.md](docs/api.md)
- **User Guide**: [docs/user-guide.md](docs/user-guide.md)
- **Developer Guide**: [docs/developer-guide.md](docs/developer-guide.md)
- **Troubleshooting**: [docs/troubleshooting.md](docs/troubleshooting.md)

## 🐛 Troubleshooting

### **Common Issues**

**503 Service Unavailable**
- Wait and retry (API may be overloaded)
- Use shorter sequences (< 400 residues)
- Switch to local ESMFold if available

**Memory Issues**
- Reduce sequence length
- Close other applications
- Use API prediction instead of local

**GPU Not Detected**
- Install CUDA-compatible PyTorch
- Verify GPU drivers
- Check VRAM requirements

### **Getting Help**
- Create an [issue](https://github.com/KeeNalytiq/ProStructure-3D/issues)
- Check [FAQ](docs/faq.md)
- Join our [Discord community](https://discord.gg/prostruct-3d)

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Meta AI** for the ESMFold model and ESM Atlas API
- **Streamlit** for the amazing web framework
- **Py3Dmol** for 3D molecular visualization
- **Biotite** for protein analysis tools
- **Plotly** for interactive visualizations

## 📚 References

- [ESM Metagenomic Atlas: The first view of the 'dark matter' of the protein universe](https://ai.facebook.com/blog/protein-folding-esmfold-metagenomics/)
- [Evolutionary-scale prediction of atomic level protein structure with a language model](https://www.biorxiv.org/content/10.1101/2022.07.20.500902v2)
- [AlphaFold's new rival? Meta AI predicts shape of 600 million proteins](https://www.nature.com/articles/d41586-022-03539-1)

## 🌟 Star History

[![Star History Chart](https://api.star-history.com/svg?repos=KeeNalytiq/ProStructure-3D&type=Date)](https://star-history.com/#KeeNalytiq/ProStructure-3D&Date)

---

<div align="center">

**Made with ❤️ by [KeeNalytiq](https://github.com/KeeNalytiq)**

[![GitHub](https://img.shields.io/badge/GitHub-100000?style=for-the-badge&logo=github&logoColor=white)](https://github.com/KeeNalytiq)
[![LinkedIn](https://img.shields.io/badge/LinkedIn-0077B5?style=for-the-badge&logo=linkedin&logoColor=white)](https://linkedin.com/in/keealytiq)
[![Twitter](https://img.shields.io/badge/Twitter-1DA1F2?style=for-the-badge&logo=twitter&logoColor=white)](https://twitter.com/keealytiq)

</div>