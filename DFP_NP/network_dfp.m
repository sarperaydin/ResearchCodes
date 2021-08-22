function w_all = network_dfp(N,own_weight,network_type)



%own_weight=self_weight;
w_com=eye(N);

if (network_type==1)
    for i=1:N
        w_com(i,mod(i,N)+1)=1;
    end
elseif(network_type==2)
     w_com(1,:)=1;
     w_com(:,1)=1;
     
else
    warning('Give 1 or 2 as input')
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
        
        if(size(inds_one,2)~=N)
            idx_own=sub2ind(size(w),inds_zero,agent*ones(1,size(inds_zero,2)));
            w(idx_own)=own_weight;
            others_weight=(1-own_weight)/(size(inds_one,2)-1);
            inds_one(inds_one==agent)=[];
        
            inds_zero_row=repmat(inds_zero,1,size(inds_one,2));
            inds_zero_cols=inds_one.*ones(size(inds_one,2),size(inds_zero,2));
            inds_zero_cols=inds_zero_cols(:)';
            idx_others=sub2ind(size(w),inds_zero_row,inds_zero_cols);
            w(idx_others)=others_weight;
        end
        w_all(:,:,agent)=sparse(w);
end


end