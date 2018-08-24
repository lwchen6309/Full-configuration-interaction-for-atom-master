function GroundState_Energy = Full_CI(Z, E_Num, basis_set)

DebugMode = false;

[Zeta,Coef,L] = ReadBasis(basis_set); 
Zeta_Num = size(Zeta,1); % Number of Gaussian basis
Orbit_Number = size(Coef,2); % Number of orbitals

% setup CSF
CSF_List = CSF_list(Orbit_Number,E_Num);
CSF_List(:,1) = [];
CSF_Number = size(CSF_List,1);

% for molecule, we should include repulsion energy of nulei
% for atom, the repulsion energy is 0, we therefore neglect it

% position of nuclei
% for atom, all position is set to [0 0 0].
% for molecule, the position depends on its geometry.
pos_atom1 = [0 0 0]; % nuclie that own wavefunction
pos_atom2 = [0 0 0]; % nuclie that interac with wavefunction
Position = repmat(pos_atom1,size(Zeta));

% One-electron integral
for ip = 1:Zeta_Num
    for iq = 1:Zeta_Num
        TmpPosition = [Position(ip,:);Position(iq,:)];
        TmpZeta = [Zeta(ip);Zeta(iq)];
        TmpL = [L(ip,:);L(iq,:)];
        % dims3 of Integral: 
        % Tab (kinetic integral)
        % Clomb (coulomb inegral)
        % Sab (overlap integral)
        one_e_integral(ip,iq,:) = One_E_Integral(TmpZeta, TmpL, TmpPosition, pos_atom2, Z);
    end
end

% Re-Normalize basis based on one-electron integral
[row,col,val] = find(Coef); % find non-zero element in coef
reNormFunc = sqrt(diag(one_e_integral(:,:,3)));
for s = 1:size(val,1)
    Coef(row(s),col(s)) = val(s) / reNormFunc(s);
end

S = Coef' * one_e_integral(:,:,3) * Coef;
% Ortho-normalize Sab , and calculate new orthcoef 
% re-calulate Tab, Coulomb on new basis (orthcoef)
Orthcoef = Coef * (S^-0.5).';
S = Orthcoef' * one_e_integral(:,:,3) * Orthcoef;
T = Orthcoef' * one_e_integral(:,:,1) * Orthcoef;
Clmb = Orthcoef' * one_e_integral(:,:,2) * Orthcoef;
h = T + Clmb; % Hamiltonian for one electron

% Two electron integral
for ip = 1:Zeta_Num
    for iq = 1:Zeta_Num
        for ir = 1:Zeta_Num
            for is = 1:Zeta_Num
                TmpPosition = [Position(ip,:);Position(iq,:);Position(ir,:);Position(is,:)];
                TmpZeta = [Zeta(ip);Zeta(iq);Zeta(ir);Zeta(is)];
                TmpL = [L(ip,:);L(iq,:);L(ir,:);L(is,:)];
                g(ip,iq,ir,is) = Two_E_Integral(TmpZeta,TmpPosition,TmpL);
            end
        end
    end
end
g = tensor(g);
g = ttm(g,{Orthcoef.',Orthcoef.',Orthcoef.',Orthcoef.'}); % Contraction by basis

% Construct Hamiltonian
for CSF1 = 1:CSF_Number
    for CSF2 = 1:CSF_Number
        Hm(CSF1,CSF2) = ...
            Hamiltonian(CSF_List(CSF1,:),CSF_List(CSF2,:));
    end
end
[~,D] = eig(Hm);
D = sort(diag(D));
GroundState_Energy = D(1);

% sub function
function H = Hamiltonian(SD1,SD2)
    % SD1 and SD2 are ON_vector (defined by second quantization)
    % Hamiltonian calculate the integral of <SD1|H|SD2> and return.
    
    if isempty(SD1)||isempty(SD2)
        return;
    end
    % Convert to CSF class : Gamma, ON_vector
    SD1 = CSF(1, SD1, Orbit_Number);
    SD2 = CSF(1, SD2, Orbit_Number);
    
    % if numbers of alpha and beta eletrons are not the same
    % return H = 0;
    if sum(SD1.Spin) ~= sum(SD2.Spin)
        H = 0;
        return;
    end
    % Difference of Ocupation Numbers
    ON_difference = sum(~ismember(SD1.ON_vector,SD2.ON_vector));
    Intersect = CSF(1,intersect(SD1.ON_vector,SD2.ON_vector),Orbit_Number);
    [I,I_Index] = setdiff(SD1.ON_vector,intersect(SD1.ON_vector,SD2.ON_vector));
    I = CSF(1,I,Orbit_Number).Orbit;
    [J,J_Index] = setdiff(SD2.ON_vector,intersect(SD1.ON_vector,SD2.ON_vector));
    J = CSF(1,J,Orbit_Number).Orbit;
    
    switch ON_difference
        case 0
            h1 = 0;
            for i =  Intersect.Orbit
                h1 = h1 + h(i,i);
            end

            % Two Electron Integral
            h2 = 0;
            for i = 1:size(Intersect.ON_vector,2)
                for j = 1:size(Intersect.ON_vector,2)
                    I = Intersect.Orbit(i);
                    J = Intersect.Orbit(j);
                    if Intersect.Spin(i) == Intersect.Spin(j)
                        h2 = h2 + g(I,I,J,J) - g(I,J,J,I);
                    else
                        h2 = h2 + g(I,I,J,J);
                    end
                end
            end

        case 1
            h1 = h(I,J) * (-1) ^ (I_Index + J_Index - 2);

            ExtendedList = [];
            h2 = 0;
            for k =  1:Orbit_Number
                [h2,ExtendedList] = Extend_DoubleExcitation(h2,SD1,SD2,[I,J,k,k],Orbit_Number,ExtendedList);
                [h2,ExtendedList] = Extend_DoubleExcitation(h2,SD1,SD2,[k,k,I,J],Orbit_Number,ExtendedList);
                [h2,ExtendedList] = Extend_DoubleExcitation(h2,SD1,SD2,[k,J,I,k],Orbit_Number,ExtendedList);
                [h2,ExtendedList] = Extend_DoubleExcitation(h2,SD1,SD2,[I,k,k,J],Orbit_Number,ExtendedList);

                SingleExcite_Integral_IKKJ = - g(I,k,k,J) * SingleExcitation(SD1,SD2,[I,J],Orbit_Number);
                h2 = h2 + SingleExcite_Integral_IKKJ;
            end

        case 2
            h1 = 0;

            ExtendedList = [];
            h2 = 0;
            [h2,ExtendedList] = Extend_DoubleExcitation(h2,SD1,SD2,[I(1),J(1),I(2),J(2)],Orbit_Number,ExtendedList);
            [h2,ExtendedList] = Extend_DoubleExcitation(h2,SD1,SD2,[I(1),J(2),I(2),J(1)],Orbit_Number,ExtendedList);
            [h2,ExtendedList] = Extend_DoubleExcitation(h2,SD1,SD2,[I(2),J(1),I(1),J(2)],Orbit_Number,ExtendedList);
            [h2,ExtendedList] = Extend_DoubleExcitation(h2,SD1,SD2,[I(2),J(2),I(1),J(1)],Orbit_Number,ExtendedList);

        otherwise
            h1 = 0;
            h2 = 0;
    end

    if DebugMode == true
        h1_Full = 0;
        for i = 1:Orbit_Number
            for j = 1:Orbit_Number
                h1_Full = h1_Full + h(i,j) * SingleExcitation(SD1,SD2,[i,j],Orbit_Number);
            end
        end

        DebugCell_Full = cell(0);
        h2_Full = 0;
        for i =  1:Orbit_Number
            for j =  1:Orbit_Number
                for k =  1:Orbit_Number
                    for l =  1:Orbit_Number
                        DoubleExcite_Integral = g(i,j,k,l) * DoubleExcitation(SD1,SD2,[i,j,k,l],Orbit_Number);
                        h2_Full = h2_Full + DoubleExcite_Integral;
                        SingleExcite_Integral = 0;
                        if j == k
                            SingleExcite_Integral = - g(i,j,k,l) * SingleExcitation(SD1,SD2,[i,l],Orbit_Number);
                            h2_Full = h2_Full + SingleExcite_Integral;
                        end
                        TmpCell = {SD1.ON_vector,SD2.ON_vector,[i,j,k,l],...
                            DoubleExcite_Integral,SingleExcite_Integral};
                        DebugCell_Full = [DebugCell_Full;TmpCell];
                    end
                end
            end
        end
        NonZeroIndex = cell2mat(DebugCell_Full(:,4)) ~= 0;
        DebugCell_Full = DebugCell_Full(NonZeroIndex,:);

        assert(abs(h1_Full - h1) < 1e-10);
        assert(abs(h2_Full - h2) < 1e-10);
    end

    H = h1 + h2/2;

    function [h2,List] = Extend_DoubleExcitation(h2,SD1,SD2,ExciationIndex,Orbit_Number,List)
        % Calculate the Double-Excitatation Integral
        % and Add to List if it's not on the List
        if isempty(List) || ~ismember(ExciationIndex,List,'Rows')
            one_e_integral = ...
                g(ExciationIndex) * DoubleExcitation(SD1,SD2,ExciationIndex,Orbit_Number);
            h2 = h2 + one_e_integral;
            List = [List; ExciationIndex];
        end
    end
end % end of Hamiltonian

end % end of Full_CI