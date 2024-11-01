# NAPDE-project
Project for the course of NAPDE of the MSc Mathematical Engineering @ PoliMi (A.Y. 2023/2024)

# Isogeometric Analysis for PDEs

## Description
This repository contains code and tools for solving partial differential equations (PDEs) using isogeometric methods implemented with the GeoPDEs library. 

## Structure
The project is organized into three main folders:
1. `geopdes`, `geopdes_hierarchical`, `nurbs-1.4.3`: those are the folders of the GeoPDEs library used for isogeometric methods. 
2. `heat_eq`: it contains the code for solving the heat equation
3. `cahn_hilliard_eq`: it contains the code for solving the Cahn-Hilliard equation, including the coupled problem with nutrients

## Installation
To clone the repository, follow these passages: 
```bash
   git clone https://github.com/reneecrispo/NAPDE-project.git

## Notes
In order to run any of the scripts, ensure that the library is set up correctly by running the file `add_path.m`