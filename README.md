# electrodes_registration

**Implementation of a High-Density EEG Cap Electrodes Registration Pipeline in MATLAB.**

## Table of Contents

- [Introduction](#introduction)
- [Installation](#installation)
- [Usage](#usage)

## Introduction

This project provides a MATLAB-based implementation of a pipeline for registering electrodes from a high-density electroencephalography (hdEEG) cap onto a subject's scalp model. Acurate electrode registration is crucial for source localization and neuroimaging analyses in neuroscience research.
 
## Installation

1. **Clone the Repository**

   ```bash
   git clone https://github.com/cernito/electrodes_registration.git
   ```
   
2. **Add to MATLAB Path**

   Open MATLAB and add the cloned directory to your path.

3. **Install Dependencies**

   Ensure that all required MATLAB toolboxes are installed

## Usage

1. **Prepare Your Data**

- Ensure you have your MRI scans, HD-EEG cap and electrode data ready.
- Supported formats:
  - MRI scan: `.nii`
  - 3D scan: `.stl`
  - Electrode data: `.elc`, `.fcsv`
 
2. **Run the Main Script**

   ```matlab
   main
   ```

3. **Follow the Pipeline**

   Firstly you will be prompted to import a MRI scan.
   The process of creation of the 3D head model takes approximately ___ second.
   After completion user is prompted to save the created file.
   
   Next you will be guided by the program to import data and interact with the process.
   You will be prompted to import a 3D head scan. 
   User will be prompted to pre-align and filter the obtained 3D scan.
   
   After completed alignment user will be prompted to import electrode data.

