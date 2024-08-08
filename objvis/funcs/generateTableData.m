function tabdat = generateTableData(condLabel)

% format: cond  r g b  balanceGroup enabled

enabled = false(size(condLabel));
enabled(2:3) = true;

r = zeros(size(condLabel));
g = r;
b = r;
r(3) = 1;
b(2) = 1;
g(4) = 0.7; b(4) = 1;
r(5) = 1; g(5) = 0.2; b(5) = 0.7;



bal = zeros(size(condLabel));
bal(2:3) = 1;
bal(4:end) = 2;


for i = 1:numel(condLabel)
    tabdat(i,:) = {condLabel{i},r(i),g(i),b(i),bal(i),enabled(i)};
end

% % 
% tabdat = {'1', 0, 0, 0, 0, false; ...
%           '2', 0, 0, 1, 1, true; ...
%           '3', 1, 0, 0, 1, true};


end