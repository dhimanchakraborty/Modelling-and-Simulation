function [t, x] = m_GillespieAlgo(react_matrix, prop_func, tspan, x0, params)

MAX_OUTPUT_LENGTH = 1000000;
num_species = size(react_matrix, 2);
T = zeros(MAX_OUTPUT_LENGTH, 1);
X = zeros(MAX_OUTPUT_LENGTH, num_species);
T(1)     = tspan(1);
X(1,:)   = x0;
rxn_count = 1;
trajectory = 3001;

while rxn_count < trajectory        
    a = prop_func(X(rxn_count,:), params);
    a0 = sum(a);
    r = rand(1,2);
    tau = -log(r(1))/a0;
    [~, mu] = histc(r(2)*a0, [0;cumsum(a(:))]);
    
    T(rxn_count+1)   = T(rxn_count)   + tau;
    X(rxn_count+1,:) = X(rxn_count,:) + react_matrix(mu,:);    
    rxn_count = rxn_count + 1;
end

t = T(1:rxn_count);
x = X(1:rxn_count,:);
disp('X');disp(x);
end