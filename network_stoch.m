function com_local = network_stoch(N,local_copies,p,eta1,eta2,actions,mode,ack_copies)

% Input: 
%local_copies, 3Local Copies
%N: number of agents
%p: success rate of communication link
%mode: if 1, limited
%      if 2, full row

com_atp=N*(N-1);

dif=actions-local_copies;
H_i=vecnorm(dif,2,2);

H_i=H_i<eta1;

com_atp=com_atp-N*(N-sum(H_i));


H_j





com_mat=zeros(N,N);
com_mat(H_i,:)=rand(sum(H_i),N);
com_mat(H_i,:)=com_mat(H_i,:)<=p;








sent_mat=com_lim(local_copies,N);

sent_local=repmat(sent_mat,1,1,N);


com_mat=repshape(com_mat,N,1,N);
com_rep=repmat(com_mat,1,N,1);


com_local=sent_local.*com_rep;



















end