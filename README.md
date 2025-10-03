# 🎈 ProStruct - 3D

Enhanced protein structure prediction app using ESMFold with robust error handling and multiple prediction methods.

[*ESMFold*](https://esmatlas.com/about) is an end-to-end single sequence protein structure predictor based on the ESM-2 language model.

## Features

- **Robust API handling** with automatic retry and fallback endpoints
- **Multiple prediction methods**: ESM Atlas API (free) or Local ESMFold (GPU required)
- **Advanced analysis**: confidence scores, hydrophobicity, charge distribution, structural metrics
- **3D visualization** with multiple coloring modes
- **Rate limiting** to prevent API overload
- **Comprehensive error handling** with helpful troubleshooting tips

## Quick Start

1. **Install dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Run the app:**
   ```bash
   streamlit run streamlit_app.py
   ```

3. **Test API connectivity (optional):**
   ```bash
   python test_api.py
   ```

## Troubleshooting 503 Errors

If you encounter "503 Service Temporarily Unavailable" errors:

### ✅ **Automatic Solutions (Built-in)**
- **Retry mechanism**: The app automatically retries with exponential backoff
- **Multiple endpoints**: Tries different API endpoints if one fails
- **Rate limiting**: Prevents overloading the API with too many requests

### 🔧 **Manual Solutions**
1. **Wait and retry**: The ESM Atlas API may be temporarily overloaded
2. **Use shorter sequences**: Sequences < 400 residues have higher success rates
3. **Switch to Local ESMFold**: If you have a GPU with 8GB+ VRAM
4. **Check internet connection**: Ensure stable connectivity

### 🖥️ **Local ESMFold Setup (Alternative)**
For offline predictions or when API is unavailable:

```bash
# Install PyTorch (GPU version recommended)
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

# Install transformers
pip install transformers

# Run with local prediction
streamlit run streamlit_app.py
# Then select "Local ESMFold (Requires GPU)" in the sidebar
```

**Requirements for Local ESMFold:**
- GPU with 8GB+ VRAM (NVIDIA recommended)
- CUDA-compatible PyTorch installation
- ~3GB disk space for model download

## API Status

The ESM Atlas API is a free service that can become overloaded during peak usage. The app includes:

- **Automatic retry logic** with exponential backoff
- **Multiple endpoint fallbacks**
- **Rate limiting** (max 1 request per 10 seconds)
- **Clear error messages** with troubleshooting tips
- **Request tracking** to monitor API usage

## Demo App

[![Streamlit App](https://img.shields.io/badge/Streamlit-FF4B4B?style=for-the-badge&logo=Streamlit&logoColor=white)](https://esmfold.streamlit.app/)

## Credit

This app was inspired by [osanseviero's app](https://huggingface.co/spaces/osanseviero/esmfold) and enhanced with robust error handling and multiple prediction methods.

## Further Reading
For more information, read the following articles:
- [ESM Metagenomic Atlas: The first view of the 'dark matter' of the protein universe](https://ai.facebook.com/blog/protein-folding-esmfold-metagenomics/)
- [Evolutionary-scale prediction of atomic level protein structure with a language model](https://www.biorxiv.org/content/10.1101/2022.07.20.500902v2)
- [AlphaFold's new rival? Meta AI predicts shape of 600 million proteins](https://www.nature.com/articles/d41586-022-03539-1)
