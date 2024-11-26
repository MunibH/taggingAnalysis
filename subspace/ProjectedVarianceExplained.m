function projve = ProjectedVarianceExplained(C, Q)
% C: Covariance matrix (n x n)
% Q: Subspace matrix (n x k), where k is the number of subspace dimensions


if isstruct(Q)
    fns = fieldnames(Q);
    for i = 1:numel(fns)
        thisfn = fns{i};
        thisQ = Q.(thisfn);
        if isstruct(C)
            thisC = C.(thisfn);
        else
            thisC = C;
        end

        C_proj = thisQ' * thisC * thisQ;
        total_variance = trace(thisC);
        projected_variance = trace(C_proj);
        projve.(thisfn) = projected_variance / total_variance;
    end
else
    C_proj = Q' * C * Q;
    total_variance = trace(C);
    projected_variance = trace(C_proj);
    projve = projected_variance / total_variance;
end
end
