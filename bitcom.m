function [com_ind, bit_inf] = bitcom(com_mat, own_ef,N,type)
%  BITCOM Provides 1-bit version of emprical frequencies
% 
% Author: Sarper Aydin
% 
% Input:
% 
% com_mat: communication matrix for a time step generated by Bernoulli distributions
% 
% own_ef: agents' own emprical freqeuncies
% 
% Output:
% 
% com_ind: Indices of agents that start to send 1-bit information
% 
% bit_inf= information package to be sent
com_ind=find(sum(com_mat,2)>0);
l_ind=length(com_ind);

if(type>-1)
    
    sent_mat=zeros(l_ind,N);
    
    for i=1:l_ind
        if(type==0)
            [~,I]=max(own_ef(i,:));
            %own_ef(i,I)=-Inf;
            sent_mat(i,I)=1;
        
        else
            Is=[];
            for j=1:type
                [M,I]=max(own_ef(i,:));
                Is=[Is I];
                own_ef(i,I)=-Inf;
                sent_mat(i,I)=M;
            end
            S=sum(sent_mat(i,:),2);
            other_val=(1-S)/(N-type);
            sent_mat(i,setdiff(1:N,Is))=other_val;
        end
    
    end
    bit_inf=sent_mat;
else
    
    bit_inf=own_ef(com_ind,:);
end


end