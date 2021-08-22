clear
clc


N=5;
%target_x_pos = [-1 1 -1 1 0]';
%target_y_pos = [-1 1 1 -1 1]';

%theta = [target_x_pos target_y_pos];

%thetas=repmat(theta,1,1,5);

%initials=rand(N,2);
%Fhat=(1/N)*rand(N,N,N);
%[si,actions,utilities] = br_target_dfp(initials,thetas,Fhat);

own_weight=0.75;
network_type=2;
w_all = network_dfp(N,own_weight,network_type);