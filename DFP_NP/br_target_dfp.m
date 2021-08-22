function [si,actions,utils] = br_target_dfp(initial_positions,target_beliefs,Fhat)


%Input:
%Initial Positions
%Targets Positions for each agent
%Local fictitious play estimates

%Output
%
N=size(initial_positions,1);
si=zeros(N,N);





%initials_rep=repmat(initial_positions,1,1,N);

%targets=reshape(target_beliefs,1,2,N);
%targets_rep=repmat(targets,N,1,1);
distance_vecs=zeros(N,2,N);
for i=1:N
    distance_vecs(:,:,i)=initial_positions(i,:)-target_beliefs(:,:,i);
end
%distance_vecs=initials_rep-targets_rep;

dists=vecnorm(distance_vecs,2,2);



%dists=permute(dists,[1 3 2]);
dists=reshape(dists,N,N,[])';

dist_inv=1./dists;

others_prob_not=zeros(N,N);
for i=1:N
    ranges=1:N;
    ranges(i)=[];
    fhat_i=Fhat(:,:,i);
    
    fhat_i=fhat_i(ranges,:);
    others_prob_not(i,:)=prod(1-fhat_i);
end


%not_selected=1-others_prob;
all_utilities=others_prob_not.*dist_inv;


%indexe gore maxlari assign et
[utils,actions]=max(all_utilities,[],2);

all_range=1:N;
actions=actions';
action_ind=sub2ind(size(si),all_range,actions);
si(action_ind)=1;

end