# electrodes_registration

**Implementation of a High-Density EEG Cap Electrodes Registration Pipeline in MATLAB.**

## Table of Contents

- [Introduction](#introduction)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Results](#results)
- [Troubleshooting](#troubleshooting)

## Introduction

This project provides a MATLAB-based implementation of a pipeline for registering electrodes from a high-density electroencephalography (hdEEG) cap onto a subject's scalp model. Acurate electrode registration is crucial for source localization and neuroimaging analyses in neuroscience research.

## Features

- **Creation of a 3D Head Model**: Creates a point cloud head model obtained from MRI scans. Created model is saved in `.stl` format.
- **Alignment of Scan With Head Model**: Aligns head scan obtained from a 3d scanner with the MRI head model.
- **Alignment of Electrodes With the HD-EEG Cap**: Aligns electrodes with 3D scan model of patient with hd-eeg cap.
- **Export Functionality**: Exports registered electrode coordinates for use in other software.

## Prerequisities 

- **MATLAB**: Version ___ or later.
- **Toolboxes**:
  - ____ 
  - ____
- **Data**"
  - MRI scans (in `.nii` format)
  - 3D scan of HD-EEG cap (in `.stl` format)
  - Electrodes positions (in `.elc` or `.fcsv` format)
 
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


## Examples

...

## Results

After running the registration pipeline, you should obtain:
- Registered electrode coordinates aligned with the patients MRI
- Visualization plots showing electrode placements.
- Exported files for use in further analyses.

## Troubleshooting

- **Common Issues**
  ...


   


