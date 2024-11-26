function proj = ProjectDataToSubspace(input_data,Q)
% input_data is (time x trials x neurons)
% Q is subspace (neurons x dims) or a struct where each field contains
% (neurons x dims) weight matrix corresponding to a subspace

if isstruct(Q)
    fns = fieldnames(Q);
    for i = 1:numel(fns)
        thisfn = fns{i};
        proj.(thisfn) = tensorprod(input_data,Q.(thisfn),3,1);
    end
else
    proj = tensorprod(input_data,Q,3,1);
end

end