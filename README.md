# GenArch

*GenArch* is pipeline written in Matlab able to compute whole genome high resolution 3D structure.



## Description

It has been conceived and developed to overcome common computational limitations typical of 3D modeling, especially in the context of single-cell Hi-C, enabling the calculation of the 3D architecture of an individual genome at a usually unattainable resolution.



## Table of Contents

- Prerequisites
- Download
- Usage
  - 1. Initial parameters
  - 2. Load necessary data
  - 3. Run the script
- Licence
- Citation




## Prerequisites

*GenArch* is a Matlab script. Therefore it requires Matlab software.
It assumes that the users has:
- Hi-C contact matrixes in the full format, at different resolutions
  - the whole genome Hi-C map at *low resolution* 
  - the intra-chromosomal Hi-C maps at *medium resolution*
  - the intra-chromosomal Hi-C maps at the resolution used to identify TAD boundaries (called *TADsRes*)
  - the Hi-C matrixes (diagonal blocks) corresponding to TADs at *high resolution*
- computed 3D structures (with any desired 3D reconstruction algorithm) in the *Nx3* matrix format, where *N* is equal to the number of points of the structure, i.e. to the number of bins in the corresponding matrix, and the three columns correspond to the 3D coordinates *x*, *y*, *z*
  - the low resolution whole-genome 3D structure
  - the medium resolution chromosomal 3D structures
  - the high resolution TADs 3D structures
- array containing the length of each chromosome in the low resolution Hi-C map (i.e. number of bins)
- array of TAD boundaries position (i.e. bins of the TADs-resolution matrix in which there is a TAD boundary), one for each chromosome




## Download

The script *GenArch.m* can easily be download from this repository.





## Usage

1. Open the script and modify the initial parameters in the first section with the values adapted for your study. 
2. Load in Matlab the necessary data.
3. Finally, run the script.

### 1. Initial parameters

In the first section there are 5 initial parameters to modify: *C*, *LowRes*, *TADsRes*, *MedRes* and *HiRes*.
Note that resolution values are expressed in *bp*, therefore *100000* indicates *100 kb*.

When opening *GenArch* for the first time, you will see the following:

```
% this section must be modified with the values adapted for the user's
% study case

C = 20; % chromosomes number
LowRes = 10000000; % enter the chosen low resolution (to compute whole-genome reference structure)
TADsRes = 2000000; % enter the chosen TADs-resolution (to identify TAD boundaries)
MedRes = 1000000; % enter the chosen medium resolution (to compute chromosomes structure)
HiRes = 100000; % enter the chosen high resolution (to compute TADs structure)

ResRatioT = TADsRes/HiRes; % scaling factor between the resolution used to identify TAD boundaries and the high resolution (used to compute TADs structure)
ResRatioM = MedRes/HiRes; % scaling factor between the medium resolution (used to compute chromosomes structure) and the high resolution (used to compute TADs structure)
```

- *C* is the number of chromosomes in the studied genome. Replace *20* with the actual number of chromosomes.
- *LowRes* is the resolution at the whole-genome level. It is the same resolution of the whole-genome Hi-C matrix (named *M_LowRes* in the script) and of the low resolution whole-genome 3D structure (called *XYZ_LowRes_gen*). Replace *10000000* with the actual chosen low resolution.
- *TADsRes* is the resolution used to identify TAD boundaries. It is the same resolution of the intra-chromosomal Hi-C matrixes (named *M_TADsRes_#c* in the script). Replace *2000000* with the actual chosen resolution.
- *MedRes* is the resolution at the chromosomal level. It is the same resolution of the chromosomal Hi-C matrixes (named *M_MedRes_#c* in the script) and of the medium resolution chromosomal 3D structures (called *XYZ_MedRes_chr#c*). Replace *1000000* with the actual chosen medium resolution.
- *HiRes* is the resolution at the TADs level. It is the same resolution of the TADs Hi-C matrixes (named *M_HiRes_#c* in the script) and of the high resolution TADs 3D structures (called *XYZ_HiRes_chr#c_TAD#t*). Replace *100000* with the actual chosen high resolution.




### 2. Load necessary data

The user must load the following data:

- *M_LowRes*: whole-genome low resolution (the one set at the *LowRes* initial parameter) Hi-C matrix, comprehensive of intra- and inter-chromosomal data, with chromosomes ordered with ascending order numbers.
- *XYZ_LowRes_gen*: low resolution (the one set at the *LowRes* initial parameter) 3D structure of the whole genome, computed with the favorite 3D reconstruction algorithm, in the *Nx3* matrix format (with *N* equal to the number of points of the structure, i.e. to the number of bins in *M_LowRes* matrix, and the three columns corresponding to the 3D coordinates *x*, *y*, *z*).
- *L_LowRes*: array of *C* (i.e. number of chromosomes) components, where each coefficient is equal to the number of bins (length) of each chromosome in the low resolution contact map *M_LowRes*.
- *M_TADsRes_#c*: intra-chromosomal contact map of each chromosome at the resolution used to compute TAD boundaries (*TADsRes*). The symbols *#c* must be replaced with the actual chromosome number (for example, replace *#c* with *1* for chromosome 1, so that *M_TADsRes_#c* becomes *M_TADsRes_1*). A total of *C* matrixes must be uploaded.
- *TB_TADsRes_#c*: array of TAD boundaries position, one for each chromosome. Each array is as long as the number of TAD boundaries identified in the corrisponding chromosome and its coefficients are equal to the number of the bins of the matrix *M_TADsRes_#c* identified as boundaries. The symbols *#c* must be replaced with the actual chromosome number (for example, replace *#c* with *1* for chromosome 1, so that *TB_TADsRes_#c* becomes *TB_TADsRes_1*). A total of *C* arrays must be uploaded.
- *M_MedRes_#c*: intra-chromosomal contact map of each chromosome at medium resolution (*MedRes*). The symbols *#c* must be replaced with the actual chromosome number (for example, replace *#c* with *1* for chromosome 1, so that *M_MedRes_#c* becomes *M_MedRes_1*). A total of *C* matrixes must be uploaded.



- *XYZ_MedRes_chr#c*: medium resolution (the one set at the *MedRes* initial parameter) 3D structure of each chromosome, computed with the favorite 3D reconstruction algorithm FROM THE MATRIX AFTER REMOVAL OF THE CENTROMERE
, in the *Nx3* matrix format (with *N* equal to the number of points of the structure, i.e. to the number of bins in *M_MedRes_#c* matrix, and the three columns corresponding to the 3D coordinates *x*, *y*, *z*). 



### 3. Run the script






## License





## Citation




