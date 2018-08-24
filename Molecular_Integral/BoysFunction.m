function fb = BoysFunction(n,x)
    %for x>0, fb: Boys function
    fb0 = @(x) erf(sqrt(x)) * sqrt(pi/4/x);
    fb = zeros(n, 1);
    if x>0
    fb(1,1) = fb0(x);
    for nn = (0:n) + 1
        if nn == n + 1
            break;
        end
        fb(nn+1,1) = (2*(nn-1) + 1)/2/x * fb(nn,1) - exp(-x)/2/x;
    end
    elseif x==0
    fb = 1./(2*(0:n).' + 1);
    end
end