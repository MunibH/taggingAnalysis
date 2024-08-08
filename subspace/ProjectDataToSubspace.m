function proj = ProjectDataToSubspace(input_data,Q)
% input_data is (time x trials x neurons)
% Q is subspace (neurons x dims)

fns = fieldnames(Q);
for i = 1:numel(fns)
    thisfn = fns{i};
    proj.(thisfn) = tensorprod(input_data,Q.(thisfn),3,1);
end

end