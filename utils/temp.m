figure('Color','w','Position',[100 100 1000 400])

%% Antidromic Collision Panel
subplot(1,2,1)
hold on; axis equal; axis off;
title('Antidromic Collision','FontSize',14)

% Soma
rectangle('Position',[0.8 1.2 0.6 0.6],'Curvature',[1 1],'FaceColor',[0.8 0.8 0.8])
text(1.1,2,'Soma','HorizontalAlignment','center')

% Axon
plot([1.4 9],[1.5 1.5],'k--','LineWidth',1)

% Stim site
plot(9,1.5,'s','MarkerSize',10,'MarkerFaceColor','r')
text(9,1.8,'Stim.','HorizontalAlignment','center')

% Spontaneous spike (soma to axon)
annotation('arrow',[0.25 0.45],[0.53 0.53],'Color','b','LineWidth',2)
text(3,1.7,'Spontaneous','Color','b')

% Antidromic spike (axon to soma)
annotation('arrow',[0.78 0.58],[0.53 0.53],'Color',[1 0.5 0],'LineWidth',2)
text(7,1.7,'Antidromic','Color',[1 0.5 0])

% Collision point
plot(5,1.5,'rx','MarkerSize',10,'LineWidth',2)
text(5,1.2,'Collision','HorizontalAlignment','center','Color','r')

%% Synaptic Transmission Panel
subplot(1,2,2)
hold on; axis equal; axis off;
title('Synaptic Transmission','FontSize',14)

% Soma
rectangle('Position',[0.8 1.2 0.6 0.6],'Curvature',[1 1],'FaceColor',[0.8 0.8 0.8])
text(1.1,2,'Soma','HorizontalAlignment','center')

% Axon
plot([1.4 9],[1.5 1.5],'k--','LineWidth',1)

% Stim site on presynaptic neuron
plot(9,3,'o','MarkerFaceColor','purple')
plot([9 9],[3 1.5],'Color','purple','LineWidth',2)
text(9,3.2,'Presynaptic','HorizontalAlignment','center','Color','purple')

% Spontaneous spike (soma to axon)
annotation('arrow',[0.25 0.45],[0.53 0.53],'Color','b','LineWidth',2)
text(3,1.7,'Spontaneous','Color','b')

% Synaptic input
annotation('arrow',[0.35 0.35],[0.72 0.58],'Color','g','LineWidth',2)
text(2.1,2,'Synaptic','Color','g')

