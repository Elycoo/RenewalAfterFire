%% plot_figures
%%
color = get(groot,'defaultAxesColorOrder');

figure(124)
hold on;

f1 = plot([0 1], [E(f_site+2) E(f_site+2)],'LineWidth',4,'Color',color(2,:));
plot([0 0 ;1 1], [E(I+2);E(I+2)],'LineWidth',4,'Color',color(1,:))
plot([0.25 0.5; 0.25 0.5], [E(f_site+2) E(f_site+2); E(I+2)],'-k','LineWidth',2)
davinci( 'arrow', 'X', [0.25 0.25], 'Y', [E(f_site+2) E(I(1)+2)], ...
    'ArrowType', 'double','Head.Width',0.05,'Shaft.Width',0.005,'Head.Length',0.2)
davinci( 'arrow', 'X', 2*[0.25 0.25], 'Y', [E(f_site+2) E(I(2)+2)] ,...
    'ArrowType', 'double','Head.Width',0.05,'Shaft.Width',0.005,'Head.Length',0.2)
text(0.75,E(1)+0.3,'\downarrow B','FontSize',24)
text(0.75,E(2)+0.3,'\downarrow E','FontSize',24)
text(0.75,E(3)+0.3,'\downarrow A','FontSize',24)
text(0.25,(E(f_site+2) + E(I(1)+2))/2,'\DeltaE_1','FontSize',24)
text(0.5 ,(E(f_site+2) + E(I(2)+2))/2,'\DeltaE_2','FontSize',24)
plot([0 1], [E(f_site+2) E(f_site+2)]+1/beta,'--k','LineWidth',2,'Color',color(5,:))
plot([0 1], [E(f_site+2) E(f_site+2)]-1/beta,'--k','LineWidth',2,'Color',color(5,:))

hold off

ylabel('"E" [a.u.]')
ylim([-1+min([E]),max([E])+1])
set(gcf,'Position', [445 253 1042 725])
set(gca,'FontSize',20)
plots = get(gca, 'Children');
l = length(plots);
legend(plots([l-2, l,l-4,1]), {'Energy level of states', 'Current state','\DeltaE between states','k_BT of system'});
% legend('Current state','Energy level of states','Delta in energies','Location','Best')
legend('Location','Best')
ax1 = gca;
ax1.XAxis.Visible = 'off';   % remove x-axis

%%
figure(123)
stairs([E E(3)],'LineWidth',4)
hold on;
% stairs([E_start E_start(3)],'--g','LineWidth',4)
plot(f_site+2.5,E(f_site+2),'.','MarkerSize',40)
plot([I+2.5; I+2.5] , [E(f_site+2) E(f_site+2); E(I+2)],'--k','LineWidth',2)
plot(I+2.5 , [E(f_site+2) E(f_site+2)],'vk','LineWidth',2)
plot(I+2.5 , E(I+2),'^k','LineWidth',2)

hold off

xticks([1.5 2.5 3.5])
xticklabels({'B','E','A'})
% ylim([-1+min([E E_start]),max([E E_start])+1])
ylim([-1+min([E]),max([E])+1])

set(gcf,'Position', [445 253 1042 725])
set(gca,'FontSize',20)

% legend('Energy level of states','Initial Energy level of states','Current state','Delta in energies','Location','Best')
legend('Energy level of states','Current state','Delta in energies','Location','Best')

%% regular field in new color
map = [[120/255   222/255   0];
       [0.9769    0.9839    0.0805];
       [0.2422    0.1504    0.6603]];
title('')
axis square
set(gca,'YDir','normal')
colormap(map)
colorbar('Ticks',[-2/3,0,2/3],...
         'TickLabels',{'B','E','A'})
set(gca,'FontSize',20)
xticks([])
yticks([])
