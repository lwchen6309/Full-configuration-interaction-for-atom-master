A FCI implementation (for atom only) in MATLAB.

(1) Setup enviroment : 
(a) MATLAB Tensor Toolbox from Sandia National Laboratory : 
http://www.sandia.gov/~tgkolda/TensorToolbox/
To run this program, please download the tensor toolbox, 
and add to path of your MATLAB.

(b) Download Basis set of atom: 
Download basis set pf atom (with format : Molpro) from EMSL Baiss Set Exchange :
link : https://bse.pnl.gov/bse/portal
paste the part "basis = { ... }" into a text file.

(2) Molecular Integral : 
A McMurchie-Davidson Algorithm is implemented to evaluate Gaussian integral.
For detail of McMurchie-Davidson Algorithm, I strongly recommend
(a) Textbook : "Trygve Helgaker et al., Molecular Electronic-Structure Theory."
(b) Short summary by Porf. Trygve Helgaker:
    (This one is important since my notation almost follow those in this slides.)
    link : http://folk.uio.no/helgaker/talks/SostrupIntegrals_10.pdf
    
(3) Execute FCI calculation :
Run Full_CI.m execute full configuration interaction and retrun ground state energy of atom.

(4) Check the calculation result
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
