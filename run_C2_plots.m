%% Plot Aircraft Weight and Thrust Required Graph

% Plot relative propoerties as a function of flight time in mins
figure
plot(FlightTime(:,1)./60,W_N(:,1)./1000) % Plot Aircraft weight
hold on
plot(FlightTime(:,1)./60,Treq_N(:,1)./1000) % Plot engine total trim thrust (requried)
hold off

% Specify figure title and axis labels
Title = 'Aircraft Weight and Engine Thrust Requirement along Flight Path'; % Title
PlotXLabel = 'Flight time elapsed, $t_{flight}$ (min)'; % SET X-AXIS LABEL
PlotYLabel = 'Aircraft Weight $W$ / Engine Thrust Required $T_{req.}$, (kN)'; % SET Y-AXIS LABEL
grid on

% Specify axis limits and ticks
xLimit(1) = 0;
xLimit(2) = 300;
yLimit(1) = 0;
yLimit(2) = 2500;
xTicks = [0 0.2 0.4 0.6 0.8 1.0].*xLimit(2);
yTicks = [0 0.2 0.4 0.6 0.8 1.0].*yLimit(2);

% Specify legend names and position
legendnames = [{'Aircraft Weight, $W$ (kN)'};...
    {'Engine Thrust, $T_{req.}$ (kN)'}]; 
legloc = 'NorthEast';

% Set fonts, line widths, colours
FONTS = 'Times New Roman';
LABSIZE = 20;
AXSIZE = 18;
LEGSIZE = 16;
LINSIZE = 2;
PAPSIZE = [14.5 14.5];
STYLES = {'-','--',':','-.','-','--',':','-.'};
COLORS = [0 0 0; 0 0 0.8; 1 0 0; 0 0.4 0; 0 0.5 1; 1 0.4 0.4; 0.5 1 0];

% Set Legend Properties
legend(gca,legendnames,'Interpreter','latex')
set(legend,'Location',legloc)
set(legend,'fontsize',LEGSIZE,'FontName',FONTS)
set(legend,'Box','off')

% Set figure properties
set(gcf,'Visible','on','PaperSize',PAPSIZE,'PaperUnits',...
    'centimeters','PaperPosition',[0.1 0.1 0.9 0.9]);
set(gcf,'Color',[1,1,1])
set(gcf,'PaperPositionMode','manual','PaperUnits','centimeters',...
    'PaperSize',PAPSIZE,'PaperPosition',[0 0 PAPSIZE(1) PAPSIZE(2)]);
set(gcf, 'Position',  [100, 100, 650, 650])

% Set axis properties
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
set(gca,'TickLength',[0.01 0.01])
set(gca,'LineWidth',1.5)
grid(gca,'off')  
box(gca,'off')
axis(gca,'square')

% Set axis labels
xlabel(PlotXLabel,'fontsize',LABSIZE,'fontname',FONTS,'Interpreter','latex')
ylabel(PlotYLabel,'fontsize',LABSIZE,'fontname',FONTS,'Interpreter','latex')
title(Title,'fontsize',LABSIZE,'fontname',FONTS,'Interpreter','latex')

% Set axis numbers
set(gca,'xlim',xLimit,'XTick',xTicks,'fontsize',AXSIZE,'fontname',FONTS)
set(gca,'ylim',yLimit,'YTick',yTicks,'fontsize',AXSIZE,'fontname',FONTS)

% Set line properties
lines = get(gca,'children');
N =  length(lines);
for b = 1:N    
    set(lines(b),'Color',COLORS(N+1-b,:),'LineWidth',LINSIZE,'LineStyle',STYLES{N+1-b});
end




