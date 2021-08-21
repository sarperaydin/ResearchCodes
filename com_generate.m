function [com_mat,com_atp] = com_generate(N,actions,local_copies,ack_copies,eta1,eta2)

atp=N*(N-1);



dif=actions-local_copies;

H_i=vecnorm(dif,2,2);









end
