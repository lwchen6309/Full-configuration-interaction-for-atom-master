function Output = Two_E_Integral(zeta,pos,L)
% two electron integral : 
% Coulomb integral

% input argument : 
% zeta, pos, L are the zeta, position and angular quantum number 
% of 4 basis, respectively.

%zeta = [a;b;c;d]
%l = [la;lb;lc;ld]
p = zeta(1) + zeta(2);
q = zeta(3) + zeta(4);
alpha = p * q / (p + q);
Pr = (zeta(1:2).' * pos(1:2,:)) ./ p;
Pq = (zeta(3:4).' * pos(3:4,:)) ./ q;
R_pq = Pr - Pq;
ab = L(1,:) + L(2,:);
cd = L(3,:) + L(4,:);
%[a;b],[A;B],[la;lb]
E1 = Eab(zeta(1:2), pos(1:2,:), L(1:2,:));
%[c;d],[C;D],[lc;ld]
E2 = Eab(zeta(3:4), pos(3:4,:), L(3:4,:));
M = R0(ab+cd, alpha, R_pq);
Output = 0;
for t1 = (0:ab(1))+1
    for u1 = (0:ab(2))+1
        for v1 = (0:ab(3))+1
            for t2 = (0:cd(1))+1
                for u2 = (0:cd(2))+1
                    for v2 = (0:cd(3))+1
                        Output = Output + ...
                            ((-1)^(t2+v2+u2-3)) * ...
                            E1(t1,u1,v1)*E2(t2,u2,v2) * ...
                            M(t1+t2-1,u1+u2-1,v1+v2-1);
                    end
                end
            end
        end
    end
end
Output = Output * 2 * (pi^2.5) / (p * q * sqrt(p + q));
end