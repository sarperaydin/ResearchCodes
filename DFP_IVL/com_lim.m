function sent_mat=com_lim(local_copies,N)



sent_mat=zeros(N,N);
for i=1:N
    [M,I]=max(local_copies(i,:,i));
    sent_mat(i,I)=M;
    other_val=(1-M)/(N-1);
    sent_mat(i,setdiff(1:N,I))=other_val;
end







end