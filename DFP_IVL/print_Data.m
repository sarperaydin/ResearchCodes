

etas_1=[0.025 0.05 0.1 0.01];
etas_2=[0.05 0.1 0.2 0.02];
etas_3=[0.01 0.02 0.04 0.004];
etas=[etas_1;etas_2;etas_3]';


%% Plotting


%% Distance to NE
figure
% plot(1:time_horizon, actions)
%%%%%%  plots
time_horizon=time_horizon_final;
%tile = tiledlayout('flow','TileSpacing','compact');
hold on
plot(1:time_horizon,dis_avg(1,:),'Color','k','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,dis_avg(2,:),'Color','r','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,dis_avg(3,:),'Color','b','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,dis_avg(4,:),'Color','g','MarkerSize',6, 'LineWidth', 2)
set(gca,'XTick',[0:(time_horizon_final/5):time_horizon],'FontSize',16 )
xlabel('Time','FontSize',16)
ax = gca;
ax.XAxis.Exponent = 2;
xticks([0:(time_horizon_final/5):time_horizon])
ylabel('Distance to equilibrium','FontSize',16)
set(gca,'FontSize',16)
grid on
set(ax, 'YScale', 'log')
set(ax,'FontSize',16)

%legend(strcat('\eta_1,\eta_2,\eta_3=',num2str(etas)),'Location','bestoutside','NumColumns',1,'FontSize',8)
legend(strcat('\eta=',num2str(etas(:,1)),'x'),'Location','southwest','NumColumns',1,'FontSize',16)
%lgd.Layout.Tile = 4;
%legend({'DFP-IVL','DFP-I','DFP-VL','DFP'},'Location','northwest','NumColumns',2)
%legend({'p=0.2-\beta=0.4','p=0.4-\beta=0.6','p=0.6-\beta=0.9'},'Location','northwest','NumColumns',1)
%legend({'DFP-IVL','DFP-I','DFP-VL','DFP'},'Location','southwest','NumColumns',2)
%legend({'1-bit','2-bit','4-bit','6-bit'},'Location','southwest','NumColumns',2)
saveas(gcf,'dist_NE.png')
print_plot('dist_nE')


%% Average Communication Error


figure
% plot(1:time_horizon, actions)
%%%%%%  plots
time_horizon=time_horizon_final;
hold on
plot(1:time_horizon,comer_avg(1,:),'Color','k','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,comer_avg(2,:),'Color','r','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,comer_avg(3,:),'Color','b','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,comer_avg(4,:),'Color','g','MarkerSize',6, 'LineWidth', 2)

set(gca,'XTick',[0:(time_horizon_final/5):time_horizon],'FontSize',16 )
xlabel('Time','FontSize',16)
ax = gca;
ax.XAxis.Exponent = 2;
xticks([0:(time_horizon_final/5):time_horizon])
ylabel('Belief Error','FontSize',16)
set(gca,'FontSize',16)
grid on
set(ax, 'YScale', 'log')
set(ax,'FontSize',16)
%legend({'p=0.2-\beta=0.3','p=0.4-\beta=0.6','p=0.6-\beta=0.9'},'Location','northwest','NumColumns',1)
%legend({'DFP-IVL','DFP-I','DFP-VL','DFP'},'Location','northwest','NumColumns',2)
saveas(gcf,'av_com.png')
print_plot('belief_er')


%% Average Acknowledgement Error

figure
% plot(1:time_horizon, actions)
%%%%%%  plots
time_horizon=time_horizon_final;
hold on
plot(1:time_horizon,acker_avg(1,:),'Color','k','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,acker_avg(2,:),'Color','r','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,acker_avg(3,:),'Color','b','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,acker_avg(4,:),'Color','g','MarkerSize',6, 'LineWidth', 2)
set(gca,'XTick',[0:(time_horizon_final/5):time_horizon],'FontSize',16 )
xlabel('Time','FontSize',16)
ax = gca;
ax.XAxis.Exponent = 2;
xticks([0:(time_horizon_final/5):time_horizon])
ylabel('Second Order Belief Error','FontSize',16)
set(gca,'FontSize',16)
grid on
set(ax, 'YScale', 'log')
set(ax,'FontSize',16)
%legend({'p=0.2-\beta=0.3','p=0.4-\beta=0.6','p=0.6-\beta=0.9'},'Location','northwest','NumColumns',1)
%legend({'DFP-IVL','DFP-I','DFP-VL','DFP'},'Location','northwest','NumColumns',2)
saveas(gcf,'ack_er.png')
print_plot('ack_er')

%% Average Communication Attepmt

figure
% plot(1:time_horizon, actions)
%%%%%%  plots
time_horizon=time_horizon_final;
hold on
plot(1:time_horizon,com_atp_avg(1,:),'Color','k','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,com_atp_avg(2,:),'Color','r','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,com_atp_avg(3,:),'Color','b','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,com_atp_avg(4,:),'Color','g','MarkerSize',6, 'LineWidth', 2)
set(gca,'XTick',[0:(time_horizon_final/5):time_horizon],'FontSize',16 )
xlabel('Time','FontSize',16)
ax = gca;
ax.XAxis.Exponent = 2;
xticks([0:(time_horizon_final/5):time_horizon])
ylabel('Communication Attempt','FontSize',16)
set(gca,'FontSize',16)
grid on
%set(ax, 'YScale', 'log')
set(ax,'FontSize',16)
%legend({'p=0.2-\beta=0.3','p=0.4-\beta=0.6','p=0.6-\beta=0.9'},'Location','northwest','NumColumns',1)
%legend({'DFP-IVL','DFP-I','DFP-VL','DFP'},'Location','northwest','NumColumns',2)
saveas(gcf,'com_atp.png')
print_plot('com_atp')