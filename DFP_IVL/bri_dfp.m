function [vec_ac,actions,utils] = bri_dfp(distances,Fhat,N)
%  BRI_DFP Provides optimal action and utility values with DFP in Target Assignment Game
% 
% Author: Sarper Aydin
% 
% Input:
% 
% Distances: distances between targets and agents positions
% 
% Fhat: local copies emprical frequencies
sz=size(Fhat);

sz3=sz(1);
others_prob_not=zeros(sz3,sz(2));
dist_inv=1./distances;
for i=1:sz3
    ranges=1:N;
    ranges(i)=[];
    fhat_i=squeeze(Fhat(i,:,ranges))';
    others_prob_not(i,:)=prod(1-fhat_i);
end
    
all_utilities=others_prob_not.*dist_inv;
[utils,actions]=max(all_utilities,[],2);
all_range=1:sz3;
actions=actions';
vec_ac=zeros(sz3,sz(2));
action_ind=sub2ind(size(vec_ac),all_range,actions);
vec_ac(action_ind)=1;
end