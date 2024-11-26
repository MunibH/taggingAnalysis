function recon = ReconstructDataFromSubspace(proj,Q)

% proj is projection of data onto subspace (time x trials x dims) or a
% struct with each field containing the same (time x trials x dims)
% Q is subspace (neurons x dims) or a struct where each field contains
% (neurons x dims) weight matrix corresponding to a subspace


if isstruct(Q)
    fns = fieldnames(Q);
    for i = 1:numel(fns)
        thisfn = fns{i};
        recon.(thisfn) = tensorprod(proj.(thisfn),Q.(thisfn),3,2);
    end
else
    recon = tensorprod(proj,Q,3,2);
end




end