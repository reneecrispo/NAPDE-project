# Isogeometric Analysis for PDEs
Project for the course of NAPDE of the MSc Mathematical Engineering @ PoliMi (A.Y. 2023/2024)

## Description
This repository contains code and tools for solving partial differential equations (PDEs) using isogeometric methods implemented with the GeoPDEs library. Our focus is on implementing and analyzing the Cahn-Hilliard equation for tumor growth. This equation will be coupled with another equation in order to model the interaction between tumor cells and nutrient concentration.

## Structure
The project is organized into three main folders:
1. `geopdes`, `geopdes_hierarchical`, `nurbs-1.4.3`: those are the folders of the GeoPDEs library used for isogeometric methods. 
2. `cahn_hilliard_equation`: it contains the code for solving the Cahn-Hilliard equation, including the coupled problem with nutrients. 

## Script
In order to run any of the scripts, ensure that the library is set up correctly by first running the file `add_path.m`. Another important scripts is `op_u_v_tp_cahn_hilliard_non_lin.m` which is the implementation of the non linear operator that can be found in the path `geopdes\inst\space\@sp_scalar`.


## Installation
To clone the repository, follow these passages: 
```bash
   git clone https://github.com/reneecrispo/NAPDE-project.git

