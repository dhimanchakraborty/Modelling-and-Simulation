function [ t, x ] = firstReactionMethod( react_matrix, propensity_fcn, tspan, x0,params)
MAX_OUTPUT_LENGTH = 1000000;
num_rxns = size(react_matrix, 1);
num_species = size(react_matrix, 2);
T = zeros(MAX_OUTPUT_LENGTH, 1);
X = zeros(MAX_OUTPUT_LENGTH, num_species);
T(1)     = tspan(1);
X(1,:)   = x0;
rxn_count = 1;
while T(rxn_count) < tspan(2)        
    a  = propensity_fcn(X(rxn_count,:), params);
    r = rand(num_rxns,1);
    taus = -log(r)./a;
    [tau, mu] = min(taus);
    disp(tau);
    T(rxn_count+1)   = T(rxn_count)   + tau;
    X(rxn_count+1,:) = X(rxn_count,:) + react_matrix(mu,:);   
    rxn_count = rxn_count + 1;         
end  
t = T(1:rxn_count);
x = X(1:rxn_count,:);disp(x);
end

