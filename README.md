# Full-configuration-interaction-for-atom
A FCI  (for atom) program implemented in MATLAB.

# Theoretical Introduction
Short introduction from wiki :
https://en.wikipedia.org/wiki/Full_configuration_interaction

(1) Full configuration interaction : 

I recommend :

  (a) a textbook introducing Quantum Chemical Calculations (mainly about electronic structure) :
  
    ref : Helgaker et al, Molecular Electronic-Structure Theory; doi:10.1002/9781119019572


  (b) a nice introduction by Gan, which shows more detail of contraction of CSF
  
    ref : The parallel implementation of a full configuration interaction program; 
  
    Gan et al, J. Chem. Phys, 119, 47 (2003); doi:10.1063/1.1575193
        
        

(2) Molecular Integral : 

A McMurchie-Davidson Algorithm is implemented to evaluate Gaussian integral.

For detail of McMurchie-Davidson Algorithm, I recommend :

  (a) Textbook : Trygve Helgaker et al., Molecular Electronic-Structure Theory.
  
  (b) Short summary by Porf. Trygve Helgaker:
  
  (This one is important since my notation almost follow those in this slides.)
      
      link : http://folk.uio.no/helgaker/talks/SostrupIntegrals_10.pdf
      
# Program
(1) Setup enviroment : 

  (a) Download MATLAB Tensor Toolbox from Sandia National Laboratory : 
  
      link : http://www.sandia.gov/~tgkolda/TensorToolbox/
  
  and add to path of your MATLAB.

  (b) Download Basis set of atom: 
  
  Download basis set pf atom (with format : Molpro) from EMSL Baiss Set Exchange :
  
      link : https://bse.pnl.gov/bse/portal
  
  paste the part "basis = { ... }" into a text file.
  
(2) Execute FCI calculation :

To FCI_demo.m contain few demo for small atom and small basis, ex : H atom / STO-3G , ... 

Full_CI.m is the main program that execute full configuration interaction and retrun ground state energy of atom.

(3) Check the calculation result

  (a) Compare the result with that from other 
  
      Quantum Chemical Calculation Software(G09, Molpro, ...)

  (b) Compare the result with online database :
  
      Computational Chemistry Comparison and Benchmark DataBase : 
      link : http://cccbdb.nist.gov/introx.asp
      go to : Calculated -> Energy -> Optimized -> Energy, and search the atom
      
Since there's no result of FCI availible, 
we compare the result with that of high level calculation : ex: CISD, CCSD(T)...
These high level calculations approch FCI nicely, 
the energy should be closed(for small atom, the energies should be the same).
