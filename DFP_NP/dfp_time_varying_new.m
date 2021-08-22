
%% DFP with time-varying network
clear
clc
rep_time=20;
N=5; % number of agents
Target_nr = N; % number of targets
time_horizon_final=50;
%time_horizon=time_horizon_final-1;
termination_time=5;
%Number of network types
total_nt=2;

%% Stats initialization
dis_avg=zeros(total_nt,time_horizon_final);
dd_avg=zeros(total_nt,time_horizon_final);
hist_cover=zeros(total_nt,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
%generate_network
%Example ring network
own_weight=0.75;
w_com=eye(N);

for i=1:N
   w_com(i,mod(i,N)+1)=1;
end

 %Assumes it is given equal weights to neighbors for each agent.
  w_all=zeros(N,N,N);
for agent = 1:N
        w=zeros(N,N);
        %w(:,agent)=ones(N,1);
       
        rang=1:N;
        inds_one=rang(w_com(agent,:)==1);
        inds_zero=rang(w_com(agent,:)~=1);
        idx_ones=sub2ind(size(w),inds_one,inds_one);
        w(idx_ones)=1;
        
        idx_own=sub2ind(size(w),inds_zero,agent*ones(1,size(inds_zero,2)));
        w(idx_own)=own_weight;
        others_weight=(1-own_weight)/(size(inds_one,2)-1);
        inds_one(inds_one==agent)=[];
        
        inds_zero_row=repmat(inds_zero,1,size(inds_one,2));
        inds_zero_cols=inds_one.*ones(size(inds_one,2),size(inds_zero,2));
        inds_zero_cols=inds_zero_cols(:)';
        idx_others=sub2ind(size(w),inds_zero_row,inds_zero_cols);
        w(idx_others)=others_weight;
        
        w_all(:,:,agent)=sparse(w);
end
%}
own_weight=0.75;
network_type=2;



%% Target assignment game
%target_x_pos = [-1 1 -1 1 0 -0.5 0.5 -0.5 0 0.5]';
%target_y_pos = [-1 1 1 -1 1 1 1 -1 -1 -1]';
%target_x_pos = [-1 1 -1 1 0]';;
%target_y_pos = [-1 1 1 -1 1]';
angles = linspace(0, 2*pi, N); % 
radius = 1;
xCenter = 0;
yCenter = 0;
target_x_pos = (radius * cos(angles) + xCenter)'; 
target_y_pos = (radius * sin(angles) + yCenter)';
%theta = [target_x_pos target_y_pos];
%theta(1,:)=0.5*theta(1,:);
%theta = theta(1:N,1:2);

belief_variance_init=0.1;
sigma_init=belief_variance_init*eye(N);
%% Agent location

%x_agent = 0.3*[-1 1 -1 1 0 -0.5 0.5 -0.5 0 0.5]';
%y_agent = 0.3*[-1 1 1 -1 1 1 1 -1 -1 -1]';
%x_agent = zeros(N,1);
%y_agent=zeros(N,1);
%action_space = 1:Target_nr;
% flag = zeros(N,1);
%% initialize belief estimates


%seed = rng; %fix random seed
%rng('default');
belief_variance=0.5;
belief_sigma=belief_variance*eye(N);
belief_mean_x=target_x_pos;
belief_mean_y=target_y_pos;






%%%%%%%%%%%%%%%%%%%%%%
%storage variables
%BTM=belief_theta_mean;
%actions = zeros(1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize histograms


 % estimated empirical frequencies stored over time
 
 
%% Start algorithm
for nt=1:total_nt

    w_all = network_dfp(N,own_weight,nt);
for rep_each=1:rep_time
    
    
    beliefs_x_init=mvnrnd(zeros(N,1),sigma_init);
    beliefs_y_init=mvnrnd(zeros(N,1),sigma_init);
    initial_positions=[beliefs_x_init', beliefs_y_init'];
    
    actions_hist=zeros(time_horizon_final,N);
    utils_hist= zeros(time_horizon_final,N);
    
    f = (1/N)*ones(N,N); % actual histogram
    Fhat=repmat(f,1,1,N);% estimated empirical frequencies
    hf= zeros(N,N,time_horizon_final);
    hFhat = zeros(N,N,time_horizon_final,N);
    beliefs_x=mvnrnd(belief_mean_x,belief_sigma,N);
    beliefs_y=mvnrnd(belief_mean_y,belief_sigma,N);
    
for t=1:time_horizon_final
    update_ratio=1/(t);
 % observe and update belief on theta
 
    if(t <= termination_time)
        beliefs_x=(1-update_ratio)*beliefs_x+update_ratio*mvnrnd(belief_mean_x,belief_sigma,N);
        beliefs_y=(1-update_ratio)*beliefs_y+update_ratio*mvnrnd(belief_mean_y,belief_sigma,N);
        target_beliefs=zeros(N,2,N);
        for i=1:N
         target_beliefs(:,:,i)=[beliefs_x(i,:)',beliefs_y(i,:)'];
        end
    end
%%% Best response
    [si,actions,utils] = br_target_dfp(initial_positions,target_beliefs,Fhat);
    
    % action and utility storage
    actions_hist(t,:)=actions;
    utils_hist(t,:)=utils;
    
    
    
    % fictitious play update
    f=(1-update_ratio)*f+update_ratio*si;
    hf(:,:,t)= f;
    
    
    for agent = 1:N
        %old_f(agent,:) = f(agent,:);
        Fhat(agent,:,agent)=f(agent,:);
    end
    
    
    
    

    %local belief update
    Fhat_2=zeros(N,N,N);
    for rec_agent=1:N
         
         for j=1:N
            f_info_j=reshape(Fhat(j,:,:),N,N)';
            Fhat_2(j,:,rec_agent)=w_all(j,:,rec_agent)*f_info_j;
         end
         
          
    end
    
    hFhat(:,:,t,:)=Fhat;
    Fhat=Fhat_2;
 kkkk=5;

end 




%% End of epoch
%%%%%%%%%%%%%%%%%%%
%% compute errors
dist = zeros(time_horizon_final,N);
for t = 1:time_horizon_final
    for agent = 1:N
        for j=1:N
            dist(t,agent) = dist(t,agent)+norm(hFhat(agent,:,t,j)-hFhat(agent,:,t,agent));
        end
    end
end
dd = sum(dist');
dd_avg(nt,:)=dd_avg(nt,:)+dd;


%% Distance to equilibrium

difs=ones(N,N)-hf(:,:,time_horizon_final);
close_NE_sol=matchpairs(difs,100);
close_NE=zeros(N,N);
idx=sub2ind(size(close_NE),close_NE_sol(:,1),close_NE_sol(:,2));
close_NE(idx)=1;

dif_NE=hf-close_NE;
dis_each=vecnorm(dif_NE,2,1);
dis=sum(dis_each,2);
dis=reshape(dis,1,time_horizon_final);
dis_avg(nt,:)=dis_avg(nt,:)+dis;

num_covered=length(unique(actions))
hist_cover(nt,num_covered)=hist_cover(nt,num_covered)+1;
end


end
%%Computing Final Stats
dis_avg=(1/rep_time)*(1/N)*dis_avg;
dd_avg=(1/(rep_time*N*(N-1)))*dd_avg;

%% FIGURES
figure
% plot(1:time_horizon, actions)
%%%%%%  plots
ccc = colormap(hsv(N));
time_horizon=time_horizon_final;
%     set(0,'DefaultAxesLineStyleOrder',{'-o'})
% for state = 1:d
hold on
for ii = 1:N
    plot(1:time_horizon,actions_hist(:,ii),'Color', ccc(ii,:),'MarkerSize',6, 'LineWidth', 2)
end

% legend('Agent 1','Agent 2', 'Agent 3', 'Agent 4', 'Agent 5', 'Average signal')
set(gca,'XTick',[0:(time_horizon_final/5):time_horizon],'FontSize',16 )
set(gca,'YTick',[1:1:N],'FontSize',16 )
%  axis([1 6 2.5 7.5])
xlabel('Time','FontSize',16)
ax = gca;
ax.XAxis.Exponent = 3\2;
xticks([0:(time_horizon_final/5):time_horizon])

ylabel('Action values','FontSize',16)
set(gca,'FontSize',16)
grid on

%%
figure
% plot(1:time_horizon, actions)
%%%%%%  plots
time_horizon=time_horizon_final;
hold on
plot(1:time_horizon,dd_avg(1,:),'Color','k','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,dd_avg(2,:),'Color','r','MarkerSize',6, 'LineWidth', 2)
set(gca,'XTick',[0:400:time_horizon],'FontSize',16 )
xlabel('Time','FontSize',16)
ax = gca;
ax.XAxis.Exponent = 2;
xticks([0:(time_horizon_final/5):time_horizon])
ylabel('Average estimation error','FontSize',16)
set(gca,'FontSize',16)
grid on
set(ax, 'YScale', 'log')
set(ax,'FontSize',16)
legend('Ring','Star')
%% 
figure
% plot(1:time_horizon, actions)
%%%%%%  plots
time_horizon=time_horizon_final;
hold on
plot(1:time_horizon,dis_avg(1,:),'Color','k','MarkerSize',6, 'LineWidth', 2)
plot(1:time_horizon,dis_avg(2,:),'Color','r','MarkerSize',6, 'LineWidth', 2)
set(gca,'XTick',[0:400:time_horizon],'FontSize',16 )
xlabel('Time','FontSize',16)
ax = gca;
ax.XAxis.Exponent = 2;
xticks([0:(time_horizon_final/5):time_horizon])
ylabel('Distance to equilibrium','FontSize',16)
set(gca,'FontSize',16)
grid on
set(ax, 'YScale', 'log')
set(ax,'FontSize',16)
legend('Ring','Star')


%% Save data

save np_dfp.mat dd_avg dis_avg hist_cover;