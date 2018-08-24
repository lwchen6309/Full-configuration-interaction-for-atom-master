function Ocu_Number = CSF_list(Orbit,E_Num)

%converter from sparse matrix to full matrix
occonv = @(x) full(sparse(1, x, 1, 1, 2*Orbit));
%create all SD
Basis = combnk(1:2*Orbit, E_Num);

GroundState_Index = 1:E_Num;
%create excitation ordered SD 
Ocu_Number = zeros(size(Basis, 1), 1);
%first coulomn of D is excitation order
for i = 1:size(Basis, 1)
    Ocu_Number(i) = sum(~ismember(GroundState_Index, Basis(i,:)));
end
Ocu_Number = [Ocu_Number, Basis];
Ocu_Number = sortrows(Ocu_Number, 1);

end