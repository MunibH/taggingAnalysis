function plotVE(ve,varargin)

if nargin > 1
    var2exp = varargin{1};
    n = numComponentsToExplainVariance(ve, var2exp);
end



f = figure;
ax = nexttile;
hold on;

xs = 1:numel(ve);
b = bar(xs,ve);

xlabel('Dim')
ylabel('%VE')

if nargin > 1
    xline(n,'k--');
    title(['# Dims for ' num2str(var2exp) ' %VE = ' num2str(n)])
end

ax.FontSize = 12;

%%

ax = nexttile;
hold on;

cve = cumsum(ve);
plot(cve,'k-','LineWidth',2);
plot(cve,'k.','MarkerSize',10)
xlabel('Dim')
ylabel('Cumulative %VE')
if nargin > 1
    xline(n,'k--');
    yline(var2exp,'k--')
end

ax.FontSize = 12;

end






