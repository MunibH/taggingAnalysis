clear,clc,close all

addpath(genpath(pwd))
%%

alm.r = 17;
alm.l = 11;
tjm1.r = 17;
tjm1.l = 8;

% Extract the data into an array
data = [alm.r, alm.l; tjm1.r, tjm1.l; sum([alm.r,tjm1.r]), sum([alm.l,tjm1.l])];

%%
close all


% Create the bar plot with stacked bars
f = figure;
f.Renderer = 'painters';
f.Position = [680   449   383   429];
ax = prettifyAxis(gca);
hBar = bar(ax, data, 'stacked');

% Set the colors
hBar(1).FaceColor = 'b'; % Blue for 'r'
hBar(2).FaceColor = 'r'; % Red for 'l'

hBar(1).EdgeColor = 'none';
hBar(2).EdgeColor = 'none';

xs = hBar(1).XData;
ys = sum(data,2);

text(xs-0.075,ys+2,num2str(ys), "FontSize",12,'FontWeight','bold')

% Set the labels and title
set(gca, 'XTickLabel', {'ALM', 'tjM1','ALM+tjM1'},'FontSize',12);
xlabel('Region','FontSize',12);
ylabel('Tagged unit count','FontSize',12);
ll = legend('Right Hemisphere', 'Left Hemisphere');
ll.Box = 'off';
ll.Location = 'northwest';
ll.FontSize = 12;




