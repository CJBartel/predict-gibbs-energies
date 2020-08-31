# predict-gibbs-energies

Repository associated with https://arxiv.org/abs/1805.08155 used to predict the temperature-dependent Gibbs energies of inorganic crystalline solids.

Also implemented in pymatgen (GibbsComputedStructureEntry) -- https://pymatgen.org/pymatgen.entries.computed_entries.html

## Data files

  ### TableS1.csv 
    
    440 compounds used for training and testing the SISSO-learned descriptor
    
  ### masses.json

    dictionary of atomic masses (amu) {element : mass}

  ### Gels.json
  
    dictionary of experimental Gibbs energies (chemical potentials) for the elements {temperature (K) : {element : G (eV/atom)}}

## Implementing descriptor

  ### PredictG.py
    
  contains class for implementing SISSO-learned descriptor; see comments within and arXiv link for details    
    
  ### POSCAR.mp-1143_Al2O3
  
  structure file for Al2O3 from Materials Project to demonstrate use of descriptor
