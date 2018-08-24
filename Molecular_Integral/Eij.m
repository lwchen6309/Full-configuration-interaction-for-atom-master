function E = Eij(zeta, pos, l, axes)
% Eij calculatethe product of 2 Gaussian, 
% and decompose it into Hermite Gaussian on 1 axis (x, y, z)

% input argument :
% zeta and pos are zeta and position of 2 Gaussian, respectively.
% l is angular quantum number, and axes is one of x,y,z coordinate(1, 2, 3, respectively)

% output argument :
% E is the amplitude of Hermite Gaussian

if l >= 0
    [a,b] = deal(zeta(1), zeta(2));
    [i,j] = deal(l(1,axes), l(2,axes));
    rab = pos(2,:) - pos(1,:);
    p = a + b;
    mu = a * b / (a + b);
    rpa = -b / p * rab;
    rpb = a / p * rab;
    Kab = exp(-mu * rab.^2);
    f(1) = rpb(axes);
    f(2) = rpa(axes);
    K = Kab(axes);
    
    % degrade of (n,t,u,v) from (i+j+k,0,0,0) to (j+k,i,0,0)
    E = zeros(i+2, j+2, i+j+2);
    E(2, 2, 2) = K;
    for hh = (0:i+j) + 2
        for r = (0:i) + 2
            for c = (0:j-1) + 2
            E(r, c+1, hh) = E(r, c, hh) * f(1) + ...
                E(r-1, c, hh) * (r-2)/2/p + ...
                E(r, c-1, hh)*(c-2)/2/p + E(r, c, hh-1)/2/p;
            end
            if r == i + 2
                break;
            end
            E(r+1, 2, hh) = E(r, 2, hh) * f(2) + ...
                E(r-1, 2, hh) * (r-2)/2/p + ...
                E(r, 1, hh) * (2-2)/2/p + E(r, 2, hh-1)/2/p;
        end
    end
    E(:,:,1) = [];
    E(1,:,:) = [];
    E(:,1,:) = [];
    E = permute(E(end,end,:),[3 1 2]);
else
    E = 0;
end
end