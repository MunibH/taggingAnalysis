function recon = ReconstructDataFromSubspace(proj,Q)

% proj is projection of data onto subspace (time x trials x dims)
% Q is subspace (neurons x dims)
recon = tensorprod(proj,Q,3,2);

end