% demo of FCI calculation on small atom
% Z is Nuclear charge, E_num is #electron of atom

% H atom
Z = 1; 
E_Num = 1;
basis_set = 'H_STO-3G.txt';
% basis_set = 'H_aTZ.txt';

% He atom
% Z = 2;
% E_Num = 2;
% basis_set = 'He_6-31G.txt';
% basis_set = 'He_STO-3G.txt';
% basis_set = 'He_DZ.txt';
% basis_set = 'He_aTZ.txt';

% Li atom
% Z = 3;
% E_Num = 3;
% basis_set = 'Li_6-31G.txt';
% basis_set = 'Li_STO-3G.txt';

E = Full_CI(Z, E_Num, basis_set);
display(E);

