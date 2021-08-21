function [tr_ind,res]=info_send(inds, com_mat,sent_mat,N)
%  INFO_SEND Brief summary of this function.
% 
% Detailed explanation of this function.
len_i=length(inds);
com_mat2=reshape(com_mat(inds,:)',N,1,len_i);
com_rep=repmat(com_mat2,1,N,1);

sent_mat2=reshape(sent_mat',1,N,len_i);
sent_rep=repmat(sent_mat2,N,1,1);
com_prod=com_rep.*sent_rep;

res=zeros(N,N,N);
res(:,:,inds)=com_prod;

tr_ind=find(res);

end