## Flexible Lightsail Simulator (2024)  
By Michael D. Kelzenberg and Ramon Gao  

### General information

  
This repository contains several versions of our code used to obtain results presented in our manuscript posted on [arXiv](https://arxiv.org/abs/2301.08894) and currently in press.

<p align="justify"> Note that two current versions of our code exist. The code in the folder "Current code (2024)" was used to obtain results for curved lightsails, while code used to study flat,
metagrating-based lightsails can be found in "Codes for manuscript".</p>

### Instructions for running simulations of flexible lightsails

<p align="justify">The main script to run for simulating flexible lightsails are files that contain "simulate" in their names. These scripts are extensively documented, and rely on calling
other scripts and functions.</p>  

<p align="justify">The first step of any simulation is to define the lightsail geometry and meshing parameters in scripts starting with "generateMesh". For example, users can define the shape 
of the lightsail (curved or flat), its geometry (radius, profile and aspect ratio for curved lightsails), the resolution of meshing (radialRings) and the type of texturing of the lightsail.
Each region/texture is associated with a specific material and type of light-matter interaction (specular reflection vs. more complex optical interaction based on look-up tables from optical
  simulations in, e.g., Lumerical or COMSOL).</p>    

<p align="justify">Material properties are defined in scripts with the name "SetupMaterialProperties". The index of a specific material is referenced in the texture sections of the mesh generation scripts.  




