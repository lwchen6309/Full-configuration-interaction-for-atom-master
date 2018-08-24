function [Zeta, Coef, L] = ReadBasis(filename)
    % read zeta, coefficient from file : 
    % zeta is a list of zeta for all basis 
    % smae as coef.
    % each basis is copied by 2*L+1 times 
    % for degeneracies of angular wavefunction
    % and L is angular quantum number.
    % ex: for   s orbital : L = [0 0 0]
    %           p orbital : L = [1 0 0], [0 1 0], [0 0 1]
    %           thus zeta, coef for p orbial is copied by 3 times, and
    %           assigned with 3 kinds of L respectively
    
    Coef = [];
    Zeta = [];
    L = [];
    if filename
        fid = fopen(filename,'r');
        % include at most s ~ h orbital
        OrbitList = {'s','p','d','f','g','h'};
        while ~feof(fid)
            ThisLine = fgetl(fid); % read the next line          
            switch ThisLine(1)
                case OrbitList;
                    [~,Index] = ismember(ThisLine(1),OrbitList);
                    AngularQuntumNum = Index - 1;
                    MaxColNum = 20;
                    ZetaFormat = ['%s %s ',repmat(' %f',1,MaxColNum)];
                    TmpZeta = textscan(ThisLine,ZetaFormat,...
                        'Delimiter',',','TreatAsEmpty','\s');
                    TmpZeta = [TmpZeta{3:end}].';
                    % Copy basis by degeneracy in Cartesian coordinate
                    L_List = Ball2Box(AngularQuntumNum,3);
                    SizeOfZeta = size(TmpZeta,1);
                    Degeneracy = size(L_List,1);
                    TmpZeta = repmat(TmpZeta,Degeneracy,1);
                    Zeta = [Zeta;TmpZeta];
                    % Copy L vector by numbers of zeta 
                    % for size conistencies
                    for i = L_List.'
                        TmpL = repmat(i.',SizeOfZeta,1);
                        L = [L;TmpL];
                    end
                case 'c'
                    % the format of Coef is [c x.x Coef]
                    % which is read as %s %f*N
                    CoefFormat = ['%s',repmat(' %f',1,MaxColNum + 1)];
                    TmpCoef = textscan(ThisLine,CoefFormat,...
                        'Delimiter',',','TreatAsEmpty','\s');
                    % TmpCoef(1:2) = {'c',n.n}
                    TmpCoef = [TmpCoef{3:end}].';
                    % Degeneracy of angular wavefunction
                    for i = 1:Degeneracy
                        Coef = blkdiag(Coef,TmpCoef);
                    end
            end
        end % end while
        fclose(fid);
        assert(size(Zeta,1) == size(L,1));
        assert(size(Zeta,1) == size(Coef,1));
    end % end if
    
end