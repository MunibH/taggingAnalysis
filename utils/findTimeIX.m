function ix = findTimeIX(t,tix, varargin)

if nargin > 2
    returnRange = varargin{1};
else
    returnRange = false;
end

for i = 1:numel(tix)
    closest_val = interp1(t,t,tix(i),'nearest');
    ix(i) = find(t==closest_val);
end

if returnRange
    ix = ix(1):ix(2);
end

end