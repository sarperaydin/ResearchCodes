clear
clc
% Time and Scale Parameters 
rep_time=1;
N=20; % number of agents
Target_nr = N; % number of targets
%rhos=[0.6 0.6 0.4 0.1]; % Update rate of EMP Histogram  5-0.8
rhos=[0.8 0.8 0.8 0.8];
p=[0.6 0.6 0.6 0.6 ]; % Communication success rate
p_ack=[0.9 0.9 0.9 0.9 ]; % Communication acknowledgement success rate
time_horizon_final=10000; % termination time 5-200
%inertias=[0.3 0.3 0.1 0.9]; %inertia paramater 5-0.1
inertias=[0.1 0.1 0.1 0.1];

%For 20 protocol comparison
%etas_1=[0.01 0.01 0 0];
%etas_2=[0.6 2 0.7 2];
%etas_3=[0.01 0.01 -1 -1];

%For 4-bit 0.002-1-0.2-rho:0.5, inertia 0.2
%

%For 20 param comparison

%etas_1=[0.0001 0.0001 0.0005 0.001];
%etas_2=[0.9 0.9 0.9 0.9];
%etas_3=[0.0001 0.0005 0.0001 0.001];
%
etas_1=[0.02 0.008 0.005 0.001];
etas_2=[0.6 0.7 0.8 0.9];
etas_3=[0.02 0.005 0.002 0.001];
etas=[etas_1;etas_2;etas_3]';
%For 20 bit comparison
%etas_1=[0.002 0.002 0.002 0.002];
%etas_2=[1 1.2 1.2 1];
%etas_3=[0.2 0.2 0.05 0.001];
%etas=[etas_1;etas_2;etas_3]';

% type indicates how many index information is used. -1 gives full
% communication.
types=[1 1 1 1];
%types=[1 1 1 1];
%Number of Different Algos for Comparison
total_nt=4;
rep=100;
%% Stats initialization
dis_avg=zeros(total_nt,time_horizon_final); % Distance to NE 
comer_avg=zeros(total_nt,time_horizon_final); % Average communication estimation Error
acker_avg=zeros(total_nt,time_horizon_final); % Average communication estimation Error
com_atp_avg=zeros(total_nt,time_horizon_final); % Average communication attempt


hist_cover=zeros(total_nt,N); % Number of covered targets



for n=1:total_nt
    
    
    for k=1:rep
    %% Target assignment game

    %the function, must be on a folder in matlab path

    %Generates position around an unit circle
        angles = 2*pi*rand(N,1); % 
        %radius = 20*rand(N,1);
        radius = randi([15 20],N,1);
        xCenter = 0;
        yCenter = 0;
        target_x_pos = (radius.*cos(angles) + xCenter);
        target_y_pos = (radius.*sin(angles) + yCenter);
        target_positions=[target_x_pos,target_y_pos];
    %%
    %Randomly generated initial positions arodun origin
        sigma_init=2*eye(N);
        x_init=mvnrnd(zeros(N,1),sigma_init)';
        y_init=mvnrnd(zeros(N,1),sigma_init)';

    %Concenated version
        initial_positions=[x_init, y_init];
        distances=zeros(N,N);
        dif=zeros(N,2);
        for i=1:N
            dif=initial_positions(i,:)-target_positions;
            distances(i,:)=vecnorm(dif,2,2);
        end

    %%%% FOR LOOP GELECEK replicate the algorithms
    % Agents' INFO kept by cell array
        agents=cell(1,4);
    %Action taken at each time step
        agents{1}=eye(N);
    % Own FP emprical frequencies
        agents{2}=(1/N)*ones(N,N);
    % Estimates on others
        agents{3}=(1/N)*ones(N,N,N);
    %agents{3}(j,:,i) gives agent j's information on i
    %Estimates on itself
        agents{4}=(1/N)*ones(N,N,N);
    %agents{4}(j,:,i) gives estimate of i on agent j's information on i
    
    % Stats per each replication
    
     hf= zeros(N,N,time_horizon_final);
       
        for t=1:time_horizon_final
            
    %% Best-response (STEP 4)
        inertia_flag=binornd(1,1-inertias(n),N,1); %% !! replace as p=1-inertia
        br_ind=find(inertia_flag);


        initial_positions(br_ind,:);
        %Fhat=reshape(agents{3}(:,:,br_ind),N-1,N,length(br_ind));
        [vec_ac,actions,utils] = bri_dfp(distances(br_ind,:),agents{3}(br_ind,:,:),N);
        agents{1}(br_ind,:)=vec_ac;
    %% Update own FP EP (STEP 5)
        
        agents{2}=(1-rhos(n))*agents{2}+rhos(n)*agents{1};
        agents{2}(agents{2}<0.0001)=0;
        hf(:,:,t)=agents{2};

%% Check conditions for communication  (STEP 6)
        %size(agents{3})
        H_ii=vecnorm(agents{2}-agents{1},2,2);

        hi_flag=H_ii>etas_1(n); % 
        hi2_flag=H_ii<etas_2(n);
        %0.01; % add smth like eta_3
        hi_flag=and(hi_flag,hi2_flag);
        
        hi_ind=find(hi_flag);
        l_hi=length(hi_ind);

        com_flag=zeros(N,N);
        for i =1:l_hi
            h_ij1=vecnorm(agents{2}(hi_ind(i),:)-agents{4}(:,:,hi_ind(i)),2,2);
            %h_ij1<eta_2;
            com_flag(hi_ind(i),:)=(h_ij1>etas_3(n))';
    
        end

        com_rate=p(n)*com_flag; %

        com_atp_avg(n,t)=com_atp_avg(n,t)+max(0,(nnz(com_rate)-N));


        com_mat=binornd(1,com_rate);
        com_mat=com_mat-diag(diag(com_mat));

        %com_ind=find(sum(com_mat,2)>0);
%% Obtain 1-bit versions of FP-EP (STEP 7)


        [com_ind, bit_inf] = bitcom(com_mat, agents{2},N,types(n));
%%
%%Commucation simulation and update (STEP 7)

%len_ci=length(com_ind);
%com_mat2=reshape(com_mat(com_ind,:)',N-1,1,len_ci);
%com_rep=repmat(com_mat2,1,N,1);
%bit_mat2=reshape(bit_inf',1,N,len_ci);
%bit_rep=repmat(bit_mat2,N-1,1,1);
%com_prod=com_rep.*bit_rep;

%agents{3}(:,:,com_ind)=com_prod;
%for i=1:len_ci
%    minus_i=1:5;
 %   minus_i(com_ind(i))=[];
 %   agents{3}(com_ind(i),:,minus_i)=reshape(com_prod(:,:,i),1,N,N-1);
%end
%size(agents{3})
%%
%%Commucation simulation and update (STEP 7)
           [tr_ind,res]=info_send(com_ind, com_mat,bit_inf,N);
           agents{3}(tr_ind)=res(tr_ind);



%% Commucation acknowledgement and update (STEP 8)

            ack_rate=com_mat*p_ack(n);
            com_ack=binornd(1,ack_rate);

            ack_ind=find(sum(com_mat,2)>0);
            len_ack=length(ack_ind);

            [tr_ind,res]=info_send(ack_ind,com_ack,bit_inf,N);
            agents{4}(tr_ind)=res(tr_ind);

            for i=1:N
                ranges=1:N;
                ranges(i)=[];
                comer_avg(n,t)=comer_avg(n,t)+sum(vecnorm(agents{2}(i,:)-agents{3}(ranges,:,i),2,2));
                acker_avg(n,t)=acker_avg(n,t)+sum(vecnorm(agents{2}(i,:)-agents{4}(ranges,:,i),2,2));
            end

 
        end


%% Distance to equilibrium

        difs=ones(N,N)-agents{2};
        close_NE_sol=matchpairs(difs,100);
        close_NE=zeros(N,N);
        idx=sub2ind(size(close_NE),close_NE_sol(:,1),close_NE_sol(:,2));
        close_NE(idx)=1;

        dif_NE=hf-close_NE;
        dis_each=vecnorm(dif_NE,2,1);
        dis=sum(dis_each,2);
        dis=reshape(dis,1,time_horizon_final);
        dis_avg(n,:)=dis_avg(n,:)+dis;

        num_covered=length(unique(actions));
        if (num_covered>0)
            hist_cover(n,num_covered)=hist_cover(n,num_covered)+1;
        end

    end
end



%% Final Averaging of Stats 

%%Computing Final Stats
dis_avg=(1/rep)*(1/N)*dis_avg;
comer_avg=(1/(rep*N*(N-1)))*comer_avg;
acker_avg=(1/(rep*N*(N-1)))*acker_avg;
com_atp_avg=(1/(rep*N*(N-1)))*com_atp_avg;


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
%legend(strcat('\eta=',num2str(etas(:,1)),'x'),'Location','best','NumColumns',1,'FontSize',16)
%lgd.Layout.Tile = 4;
%legend({'DFP-IVL','DFP-I','DFP-VL','DFP'},'Location','northwest','NumColumns',2)
%legend({'p=0.2-\beta=0.4','p=0.4-\beta=0.6','p=0.6-\beta=0.9'},'Location','northwest','NumColumns',1)
%legend({'DFP-IVL','DFP-VL','DFP-I','DFP'},'Location','southwest','NumColumns',2)
name1=strcat(num2str(etas_1(1)),'<h_{ii}<',num2str(etas_2(1)),'; h_{ij}>',num2str(etas_3(1)));
name2=strcat(num2str(etas_1(2)),'<h_{ii}<',num2str(etas_2(2)),'; h_{ij}>',num2str(etas_3(2)));
name3=strcat(num2str(etas_1(3)),'<h_{ii}<',num2str(etas_2(3)),'; h_{ij}>',num2str(etas_3(3)));
name4=strcat(num2str(etas_1(4)),'<h_{ii}<',num2str(etas_2(4)),'; h_{ij}>',num2str(etas_3(4)));
legend({name1,name2,name3,name4},'Location','southwest','NumColumns',1)
%legend({'4-bit','6-bit'},'Location','southwest','NumColumns',2)
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


%%Save data

FileName1=['Dist_ne-','Agents',num2str(N),'-',datestr(now,'mmmmddyyyy-HHMMAM'),'.mat'];
FileName2=['Com_er-','Agents',num2str(N),'-',datestr(now,'mmmmddyyyy-HHMMAM'),'.mat'];
FileName3=['Ack_er-','Agents',num2str(N),'-',datestr(now,'mmmmddyyyy-HHMMAM'),'.mat'];
FileName4=['Cover','Agents',num2str(N),'-',datestr(now,'mmmmddyyyy-HHMMAM'),'.mat'];
FileName5=['Com_atp','Agents',num2str(N),'-',datestr(now,'mmmmddyyyy-HHMMAM'),'.mat'];

save(FileName1,'dis_avg');
save(FileName2,'comer_avg');
save(FileName3,'acker_avg');
save(FileName4,'hist_cover');
save(FileName5,'com_atp_avg');

%save dis_avg;
%save comer_avg;
%save acker_avg;
%save hist_cover;









