function [dydt] = conditionalODE(t,y,...
    coalescent_rates,sum,connectivity)
% does the simulation of joint line state probs
    
    dydt = zeros(size(sum));    
    
    % state change from coalescing
    for i = 1 : length(dydt)
        dydt(i) = -sum(i)*coalescent_rates*y(i);
    end
    
    % state change from migration
    dydt = dydt + (y'*connectivity)';
end