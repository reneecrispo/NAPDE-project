# Isogeometric Analysis for PDEs
Project for the course of NAPDE of the MSc Mathematical Engineering @ PoliMi (A.Y. 2023/2024)

## Description
This repository contains code and tools for solving partial differential equations (PDEs) using isogeometric methods implemented with the GeoPDEs library. Our focus is on implementing and analyzing the Cahn-Hilliard equation for tumor growth. This equation will be coupled with another equation in order to model the interaction between tumor cells and nutrient concentration.

## Structure
The project is organized into three main folders:
1. `geopdes`, `geopdes_hierarchical`, `nurbs-1.4.3`: those are the folders of the GeoPDEs library used for isogeometric methods. 
2. `heat_equation`: it contains the code for solving the heat equation
3. `cahn_hilliard_equation`: it contains the code for solving the Cahn-Hilliard equation, including the coupled problem with nutrients. It also contains a subfolder `results` where all the VTS files with simulation results are stored for each time instance.

## Script
In order to run any of the scripts, ensure that the library is set up correctly by first running the file `add_path.m`.

## Installation
To clone the repository, follow these passages: 
```bash
   git clone https://github.com/reneecrispo/NAPDE-project.git

