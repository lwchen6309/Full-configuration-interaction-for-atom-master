classdef CSF
    properties
        Amplitude; % amplitude of CSF (Configuration State Function)
        ON_vector; % occupation number (defined by second quantization)
        OrbitNumber; % total number of orbits in the basis
        Orbit; % index of orbit
        Spin; % spin of electron : beta spin is 1 and alpha spin is 0
    end
    methods
    % The method define the following function : 
    % Construction of CSF object
    % The Dot product of CSF object
    % Creation and Annihilation Operator
    % Spin-specified Excitation Operator
    % Spin-unspecified Single and Double Excitation Operator
    
    % Constructor
    function obj = CSF(Amplitude,ON_vector,OrbitNumber)
        if nargin > 0
            assert(logical(prod(ON_vector > 0)));
            obj = CSF;
            obj.Amplitude = Amplitude;
            obj.ON_vector = ON_vector;
            obj.OrbitNumber = OrbitNumber;
            obj.Spin = obj.ON_vector > OrbitNumber;%beta = 1, alpha = 0
            obj.Orbit = obj.ON_vector - obj.Spin * OrbitNumber;
        end
    end
    
    % Display Amplitude and ON_vector
    function show(CSF)
        display([int2str(CSF.Amplitude),' | ',int2str(CSF.ON_vector),' >']);
    end
    
    function u = Creation(q,v)
        % The input argument q is a spin-orbit (an integer),...
        % and v is a CSF
        % This function add q in v if q is not in v,
        % and return the resulted CSF u
        if isempty(q) || v.Amplitude == 0
            u = v;
        else
            if ~ismember(q,v.ON_vector)
                ON_vector = sort([v.ON_vector,q]);
                [~,Index_of_p] = ismember(q,ON_vector);
                Amplitude = v.Amplitude * (-1) ^ (Index_of_p - 1);
                u = CSF(Amplitude,ON_vector,v.OrbitNumber);
            else
                u = v;
                u.Amplitude = 0;
            end
        end
    end

    function u = Annihilation(p,v)
        % The input argument q is a spin-orbit (an integer),...
        % and v is a CSF (an integer array)
        % This function remove q in v if q is in v,
        % and return the resulted CSF u
        if isempty(p) || v.Amplitude == 0
            u = v;
        else
            if ismember(p,v.ON_vector)
                [~,Index_of_p] = ismember(p,v.ON_vector);
                ON_vector = v.ON_vector;
                ON_vector(Index_of_p) = [];
                Amplitude = v.Amplitude * (-1) ^ (Index_of_p - 1);
                u = CSF(Amplitude,ON_vector,v.OrbitNumber);
            else
                u = v;
                u.Amplitude = 0;
            end
        end
    end
    
    % Spin-specified Excitation Operator
    function u = E(i,j,v)
        u = Annihilation(j,v);
        u = Creation(i,u);
    end % end of E
    
    % Dot product of two CSF overloading to mtimes (Multiply)
    function DotProduct = mtimes(ON1,ON2)
        if size(ON1.ON_vector) == size(ON2.ON_vector)
            if ON1.ON_vector == ON2.ON_vector
                DotProduct = 1;
            else
                DotProduct = 0;
            end
            DotProduct = DotProduct * ON1.Amplitude * ON2.Amplitude;
        else
            DotProduct = 0;
        end
    end
    
    % Spin-unspecified Single and Double Excitation Operator
    function Integral = SingleExcitation(SD1,SD2,ExciationIndex,Orbit_Number)
        assert(size(ExciationIndex,2) == 2);
        [i,j] = deal(ExciationIndex(1),ExciationIndex(2));
        Integral = SD1 * E(i,j,SD2) + ...
            SD1 * E(i + Orbit_Number,j + Orbit_Number,SD2);
    end
    
    function Integral = DoubleExcitation(SD1,SD2,ExciationIndex,Orbit_Number)
        assert(size(ExciationIndex,2) == 4);
        [i,j,k,l] = deal(ExciationIndex(1),ExciationIndex(2),ExciationIndex(3),ExciationIndex(4));
        Integral = SD1 * E(i,j,E(k,l,SD2)) +...
            SD1 * E(i,j,E(k + Orbit_Number,l + Orbit_Number,SD2)) +...
            SD1 * E(i + Orbit_Number,j + Orbit_Number,E(k,l,SD2)) +...
            SD1 * E(i + Orbit_Number,j + Orbit_Number,E(k + Orbit_Number,l + Orbit_Number,SD2));
    end
    
    end % end of methods
end % end of class