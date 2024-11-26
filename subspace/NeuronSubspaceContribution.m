function neuron_contrib = NeuronSubspaceContribution(Q,FR)
% This function calculates the contribution of each neuron to the subspace
% Inputs:
%   FR   - firing rate of neurons (neurons,1)
%   Q    - matrix defining the subspace, of size (neurons, dimensions)



if isstruct(Q)
    fns = fieldnames(Q);
    for i = 1:numel(fns)
        thisfn = fns{i};
        thisQ = Q.(thisfn);
        
        W = sqrt(sum(thisQ.^2,2));

        % c = W .* FR;
        c = W;
        neuron_contrib.(thisfn) = c / sum(c);
    end
else
    W = sqrt(sum(Q.^2,2));

    % c = W .* FR;
    c = W;
    neuron_contrib = c / sum(c);
end


end
