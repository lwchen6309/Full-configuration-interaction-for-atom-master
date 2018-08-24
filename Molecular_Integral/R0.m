function Rn = R0(l,p,R)
% Hermite Coulomb integral

% input argument :
% l is angular momentum
% p is total exponent of 4 basis
% R is sub-index = [t, u, v]

% Construct Rntuv in v-axes
% build R(n,0,0,0) from n = 0 to i+j+k
Rn0 = (-2*p).^(0:sum(l,2)).' .* BoysFunction(sum(l,2),p*norm(R)^2);
% degrade of (n, t, u, v) from (0~i+j+k,0,0,0) to (0~j+k,0,0,0~k)
Rnv = zeros(l(3)+2,1,sum(l,2)+1);
Rnv(2,1,:) = Rn0;
for hh = (0:l(3)-1)+2
    temp = R(3)*Rnv(hh,1,:) + (hh-2)*Rnv(hh-1,1,:);
    temp = cat(3,temp(:,:,2:end),0);
    Rnv(hh+1,1,:) = temp;
end
Rnv(1,:,:) = [];
Rnv = Rnv(:,:,1:l(1)+l(2)+1);
Rnv = permute(Rnv,[3 1 2]);
Rn = [];
for ii = Rnv
    Rn = cat(3,Rn,degrade(ii));
end

% Construct Rntuc in t,u axes
function Rn = degrade(Rn0)
    % degrade of (n,t,u,v) from (i+j+k,0,0,0) to (j+k,i,0,0)
    i1 = l(1);
    j1 = l(2);
    
    Rn = zeros(i1+2,j1+2,i1+j1+1);
    Rn(2,2,:) = Rn0;
    for r = (0:i1-1)+2
        % [a b c] => [b c 0] in n-axes
        tmpl = R(1)*Rn(r,2,:) + (r-2)*Rn(r-1,2,:);
        tmpl = cat(3,tmpl(:,:,2:end),0);
        Rn(r+1,2,:) = tmpl;
    end
    Rn(1,:,:) = [];
    Rn = Rn(:,:,1:j1+1);
    for c = (0:j1-1)+2
        tmpl = R(2)*Rn(:,c,:) + (c-2)*Rn(:,c-1,:);
        tmpl = cat(3,tmpl(:,:,2:end),zeros(size(tmpl,1),1));
        Rn(:,c+1,:) = tmpl;
    end
    Rn(:,1,:) = [];
    Rn = Rn(:,:,1);
end
end