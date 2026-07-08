# DenseArrayToolkit-GUI
[![License](https://img.shields.io/badge/License-GPL--3.0-blue)](https://www.gnu.org/licenses/gpl-3.0.en.html#license-text)
![Last Commit](https://img.shields.io/github/last-commit/seismic-chen/DenseArrayToolkit-GUI)
![GitHub Pull Requests](https://img.shields.io/github/issues-pr/seismic-chen/DenseArrayToolkit-GUI)
![GitHub Issues](https://img.shields.io/github/issues/seismic-chen/DenseArrayToolkit-GUI)
![Tests](https://img.shields.io/badge/Tests-passing-brightgreen)
![Code Coverage](https://img.shields.io/badge/coverage-30%25-yellow)
![Code Size](https://img.shields.io/github/languages/code-size/seismic-chen/DenseArrayToolkit-GUI)
![Repo Size](https://img.shields.io/github/repo-size/seismic-chen/DenseArrayToolkit-GUI)
![Docs](https://img.shields.io/badge/DAT%20docs-undergoing-brightgreen)
[![DOI](https://img.shields.io/badge/DOI-10.5281/zenodo.21186268-pink)](https://doi.org/10.5281/zenodo.21186268)

![GitHub Stars](https://img.shields.io/github/stars/seismic-chen/DenseArrayToolkit-GUI?style=social)
![GitHub Forks](https://img.shields.io/github/forks/seismic-chen/DenseArrayToolkit-GUI?style=social)


## DAT-GUI

**DAT-GUI** is an open-source MATLAB package for receiver function processing and imaging of dense seismic array data.

## Requirements

- [MATLAB](https://www.mathworks.com/products/matlab.html) **R2021a** or later

## Installation

### 1. Get the software

Clone the repository:

```bash
git clone https://github.com/seismic-chen/DenseArrayToolkit-GUI.git
```

Or, download the ZIP file from GitHub and unzip it.

### 2. Download example data

Download the example dense-array datasets from Zenodo:

- **DOI:** [10.5281/zenodo.21186268](https://doi.org/10.5281/zenodo.21186268)

Two test datasets are included:

| Folder | Array type |
|--------|------------|
| `DAT_event_waveforms_BY` | Distributed array |
| `DAT_event_waveforms_SL` | Linear array |

## Quick Start

1. Open MATLAB.

2. Change to the main directory of the repository:

   ```matlab
   cd /path/to/DenseArrayToolkit-GUI
   ```

3. Start the GUI from the MATLAB Command Window:

   ```matlab
   DenseArrayToolkit_GUI
   ```

4. Load the example data and try the processing workflow.

## Modifying the Source Code

To modify the source code or add new features:

1. Open `*.mlapp` in MATLAB App Designer.
2. Switch to **Code View** to edit the code.

## Overview and main interface

<p align="center">
  <img src="https://github.com/user-attachments/assets/f4e895ae-4531-4daa-a364-d544beebfeeb" alt="DAT-GUI main interface" width="613">
</p>
<p align="center"><em>Figure 1. Processing workflow in DAT-GUI.</em></p>

<p align="center">
  <img src="https://github.com/user-attachments/assets/de851fb2-a1d5-4155-9e58-ce95dff4d466" alt="DAT-GUI processing workflow" width="613">
</p>
<p align="center"><em>Figure 2. DAT-GUI main interface.</em></p>




