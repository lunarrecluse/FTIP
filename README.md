# FTIP
A method for global protein surface comparison called FTIP, Furthest point sampling (FPS)-enhanced Triangulation-based Iterative closest point (ICP) for Protein surface comparison (PSC), which is applied to classifying proteins using only their surface shape information.
This package provides an implementation of the inference pipeline of FTIP.
![Picture1](https://github.com/lunarrecluse/FTIP/assets/9223804/50babc33-645f-4b9b-a9a2-bbb6d66f584d)

Please refer to the FTIP paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7214043/) for a detailed method description.
1. Data Preparation
   Input files must be in the PDB format. It's better to use NACCESS to 
   generate surface atoms and all the atoms with SASA=0.0 should be removed.
   Some example input files are in test_surface/

2. Running extraction script
   (1) FPS
   python extract_atom_FPS.py [Input PDB-surfae] [the number of feature points] [mode]

   Description:
   Input PDB-surfae: the surface atoms which in PDB format
   the number of feature points: The number of output feature points
   mode: use MAX for FPS_max, and MIN for FPS_min

   Example:
   python extract_atom_FPS.py 1a2pC-surf.asa 30 MIN
   python extract_atom_FPS.py 1a2pC-surf.asa 30 MAX

   (2) Kmean
   python extract_atom_kmean.py [Input PDB-surfae] [the number of feature points]

   Example:
   python extract_atom_kmean.py 1a2pC-surf.asa 30

3. Running two-steps extraction script
   FPS
   python extract_atom_FPS_EM.py [Input PDB-surfae] [the number of feature points] [grid size] [mode]
   
   Description:
   Input PDB-surfae: the surface coordinate which in x3d format
   the number of feature points: The number of output feature points
   grid size: first step grid size
   mode: use MAX for FPS_max, and MIN for FPS_min


4. Run surface comparison 
   Edit config.txt in build folder in protein_registration_version1.2/protein_registration/build/
   Under build/ folder:
   ./protein_registration
