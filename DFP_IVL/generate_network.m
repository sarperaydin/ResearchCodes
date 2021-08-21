D=[1,1,1,1,1];
%% Time varying circle network
w1(:,:,1)=[1,0,0,0,0,0,0,0,0,0;
    0.25,0.75,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w1(:,:,2)= [1,0,0,0,0,0,0,0,0,0;
    0,0.75,0.25,0,0,0,0,0,0,0;
    0,0.25,0.75,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w1(:,:,3)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,0.75,0.25,0,0,0,0,0,0;
    0,0,0.25,0.75,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w1(:,:,4)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,0.75,0.25,0,0,0,0,0; 
    0,0,0,0.25,0.75,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w1(:,:,5)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,0.75,0.25,0,0,0,0;
    0,0,0,0,0.25,0.75,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w1(:,:,6)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,0.75,0.25,0,0,0;
    0,0,0,0,0,0.25,0.75,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w1(:,:,7)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,0.75,0.25,0,0;
    0,0,0,0,0,0,0.25,0.75,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w1(:,:,8)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,0.75,0.25,0;
    0,0,0,0,0,0,0,0.25,0.75,0;
    0,0,0,0,0,0,0,0,0,1];
w1(:,:,9)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,0.75,0.25;
    0,0,0,0,0,0,0,0,0.25,0.75];
w1(:,:,10)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0.25,0,0,0,0,0,0,0,0,0.75];
for i=1:time_horizon_final
   w1(:,:,i+10)=w1(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w2(:,:,1)=[0.75,0.25,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w2(:,:,2)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0.25,0.75,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w2(:,:,3)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,0.75,0.25,0,0,0,0,0,0;
    0,0,0.25,0.75,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w2(:,:,4)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,0.75,0.25,0,0,0,0,0; 
    0,0,0,0.25,0.75,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w2(:,:,5)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,0.75,0.25,0,0,0,0;
    0,0,0,0,0.25,0.75,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w2(:,:,6)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,0.75,0.25,0,0,0;
    0,0,0,0,0,0.25,0.75,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w2(:,:,7)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,0.75,0.25,0,0;
    0,0,0,0,0,0,0.25,0.75,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w2(:,:,8)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,0.75,0.25,0;
    0,0,0,0,0,0,0,0.25,0.75,0;
    0,0,0,0,0,0,0,0,0,1];
w2(:,:,9)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,0.75,0.25;
    0,0,0,0,0,0,0,0,0.25,0.75];
w2(:,:,10)= [0.75,0,0,0,0,0,0,0,0,0.25;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0.25,0,0,0,0,0,0,0,0,0.75];
for i=1:time_horizon_final
   w2(:,:,i+10)=w2(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w3(:,:,1)=[0.75,0.25,0,0,0,0,0,0,0,0;
    0.25,0.75,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w3(:,:,2)= [1,0,0,0,0,0,0,0,0,0;
    0,0.75,0.25,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w3(:,:,3)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0.25,0.75,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w3(:,:,4)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,0.75,0.25,0,0,0,0,0; 
    0,0,0,0.25,0.75,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w3(:,:,5)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,0.75,0.25,0,0,0,0;
    0,0,0,0,0.25,0.75,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w3(:,:,6)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,0.75,0.25,0,0,0;
    0,0,0,0,0,0.25,0.75,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w3(:,:,7)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,0.75,0.25,0,0;
    0,0,0,0,0,0,0.25,0.75,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w3(:,:,8)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,0.75,0.25,0;
    0,0,0,0,0,0,0,0.25,0.75,0;
    0,0,0,0,0,0,0,0,0,1];
w3(:,:,9)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,0.75,0.25;
    0,0,0,0,0,0,0,0,0.25,0.75];
w3(:,:,10)= [0.75,0,0,0,0,0,0,0,0,0.25;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0.25,0,0,0,0,0,0,0,0,0.75];
for i=1:time_horizon_final
   w3(:,:,i+10)=w3(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w4(:,:,1)=[0.75,0.25,0,0,0,0,0,0,0,0;
    0.25,0.75,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w4(:,:,2)= [1,0,0,0,0,0,0,0,0,0;
    0,0.75,0.25,0,0,0,0,0,0,0;
    0,0.25,0.75,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w4(:,:,3)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,0.75,0.25,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w4(:,:,4)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0.25,0.75,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w4(:,:,5)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,0.75,0.25,0,0,0,0;
    0,0,0,0,0.25,0.75,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w4(:,:,6)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,0.75,0.25,0,0,0;
    0,0,0,0,0,0.25,0.75,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w4(:,:,7)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,0.75,0.25,0,0;
    0,0,0,0,0,0,0.25,0.75,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w4(:,:,8)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,0.75,0.25,0;
    0,0,0,0,0,0,0,0.25,0.75,0;
    0,0,0,0,0,0,0,0,0,1];
w4(:,:,9)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,0.75,0.25;
    0,0,0,0,0,0,0,0,0.25,0.75];
w4(:,:,10)= [0.75,0,0,0,0,0,0,0,0,0.25;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0.25,0,0,0,0,0,0,0,0,0.75];
for i=1:time_horizon_final
   w4(:,:,i+10)=w4(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w5(:,:,1)=[0.75,0.25,0,0,0,0,0,0,0,0;
    0.25,0.75,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w5(:,:,2)= [1,0,0,0,0,0,0,0,0,0;
    0,0.75,0.25,0,0,0,0,0,0,0;
    0,0.25,0.75,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w5(:,:,3)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,0.75,0.25,0,0,0,0,0,0;
    0,0,0.25,0.75,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w5(:,:,4)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,0.75,0.25,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w5(:,:,5)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0.25,0.75,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w5(:,:,6)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,0.75,0.25,0,0,0;
    0,0,0,0,0,0.25,0.75,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w5(:,:,7)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,0.75,0.25,0,0;
    0,0,0,0,0,0,0.25,0.75,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w5(:,:,8)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,0.75,0.25,0;
    0,0,0,0,0,0,0,0.25,0.75,0;
    0,0,0,0,0,0,0,0,0,1];
w5(:,:,9)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,0.75,0.25;
    0,0,0,0,0,0,0,0,0.25,0.75];
w5(:,:,10)= [0.75,0,0,0,0,0,0,0,0,0.25;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0.25,0,0,0,0,0,0,0,0,0.75];
for i=1:time_horizon_final
   w5(:,:,i+10)=w5(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w6(:,:,1)=[0.75,0.25,0,0,0,0,0,0,0,0;
    0.25,0.75,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w6(:,:,2)= [1,0,0,0,0,0,0,0,0,0;
    0,0.75,0.25,0,0,0,0,0,0,0;
    0,0.25,0.75,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w6(:,:,3)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,0.75,0.25,0,0,0,0,0,0;
    0,0,0.25,0.75,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w6(:,:,4)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,0.75,0.25,0,0,0,0,0; 
    0,0,0,0.25,0.75,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w6(:,:,5)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0.25,0.75,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w6(:,:,6)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0.25,0.75,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w6(:,:,7)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,0.75,0.25,0,0;
    0,0,0,0,0,0,0.25,0.75,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w6(:,:,8)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,0.75,0.25,0;
    0,0,0,0,0,0,0,0.25,0.75,0;
    0,0,0,0,0,0,0,0,0,1];
w6(:,:,9)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,0.75,0.25;
    0,0,0,0,0,0,0,0,0.25,0.75];
w6(:,:,10)= [0.75,0,0,0,0,0,0,0,0,0.25;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0.25,0,0,0,0,0,0,0,0,0.75];
for i=1:time_horizon_final
   w6(:,:,i+10)=w6(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w7(:,:,1)=[0.75,0.25,0,0,0,0,0,0,0,0;
    0.25,0.75,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w7(:,:,2)= [1,0,0,0,0,0,0,0,0,0;
    0,0.75,0.25,0,0,0,0,0,0,0;
    0,0.25,0.75,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w7(:,:,3)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,0.75,0.25,0,0,0,0,0,0;
    0,0,0.25,0.75,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w7(:,:,4)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,0.75,0.25,0,0,0,0,0; 
    0,0,0,0.25,0.75,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w7(:,:,5)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,0.75,0.25,0,0,0,0;
    0,0,0,0,0.25,0.75,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w7(:,:,6)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,0.75,0.25,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w7(:,:,7)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0.25,0.75,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w7(:,:,8)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,0.75,0.25,0;
    0,0,0,0,0,0,0,0.25,0.75,0;
    0,0,0,0,0,0,0,0,0,1];
w7(:,:,9)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,0.75,0.25;
    0,0,0,0,0,0,0,0,0.25,0.75];
w7(:,:,10)= [0.75,0,0,0,0,0,0,0,0,0.25;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0.25,0,0,0,0,0,0,0,0,0.75];
for i=1:time_horizon_final
   w7(:,:,i+10)=w7(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w8(:,:,1)=[0.75,0.25,0,0,0,0,0,0,0,0;
    0.25,0.75,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w8(:,:,2)= [1,0,0,0,0,0,0,0,0,0;
    0,0.75,0.25,0,0,0,0,0,0,0;
    0,0.25,0.75,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w8(:,:,3)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,0.75,0.25,0,0,0,0,0,0;
    0,0,0.25,0.75,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w8(:,:,4)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,0.75,0.25,0,0,0,0,0; 
    0,0,0,0.25,0.75,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w8(:,:,5)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,0.75,0.25,0,0,0,0;
    0,0,0,0,0.25,0.75,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w8(:,:,6)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,0.75,0.25,0,0,0;
    0,0,0,0,0,0.25,0.75,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w8(:,:,7)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,0.75,0.25,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w8(:,:,8)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0.25,0.75,0;
    0,0,0,0,0,0,0,0,0,1];
w8(:,:,9)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,0.75,0.25;
    0,0,0,0,0,0,0,0,0.25,0.75];
w8(:,:,10)= [0.75,0,0,0,0,0,0,0,0,0.25;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0.25,0,0,0,0,0,0,0,0,0.75];
for i=1:time_horizon_final
   w8(:,:,i+10)=w8(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w9(:,:,1)=[0.75,0.25,0,0,0,0,0,0,0,0;
    0.25,0.75,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w9(:,:,2)= [1,0,0,0,0,0,0,0,0,0;
    0,0.75,0.25,0,0,0,0,0,0,0;
    0,0.25,0.75,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w9(:,:,3)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,0.75,0.25,0,0,0,0,0,0;
    0,0,0.25,0.75,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w9(:,:,4)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,0.75,0.25,0,0,0,0,0; 
    0,0,0,0.25,0.75,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w9(:,:,5)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,0.75,0.25,0,0,0,0;
    0,0,0,0,0.25,0.75,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w9(:,:,6)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,0.75,0.25,0,0,0;
    0,0,0,0,0,0.25,0.75,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w9(:,:,7)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,0.75,0.25,0,0;
    0,0,0,0,0,0,0.25,0.75,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w9(:,:,8)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,0.75,0.25,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w9(:,:,9)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0.25,0.75];
w9(:,:,10)= [0.75,0,0,0,0,0,0,0,0,0.25;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0.25,0,0,0,0,0,0,0,0,0.75];
for i=1:time_horizon_final
   w9(:,:,i+10)=w9(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w10(:,:,1)=[0.75,0.25,0,0,0,0,0,0,0,0;
    0.25,0.75,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w10(:,:,2)= [1,0,0,0,0,0,0,0,0,0;
    0,0.75,0.25,0,0,0,0,0,0,0;
    0,0.25,0.75,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w10(:,:,3)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,0.75,0.25,0,0,0,0,0,0;
    0,0,0.25,0.75,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w10(:,:,4)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,0.75,0.25,0,0,0,0,0; 
    0,0,0,0.25,0.75,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w10(:,:,5)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,0.75,0.25,0,0,0,0;
    0,0,0,0,0.25,0.75,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w10(:,:,6)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,0.75,0.25,0,0,0;
    0,0,0,0,0,0.25,0.75,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w10(:,:,7)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,0.75,0.25,0,0;
    0,0,0,0,0,0,0.25,0.75,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
w10(:,:,8)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,0.75,0.25,0;
    0,0,0,0,0,0,0,0.25,0.75,0;
    0,0,0,0,0,0,0,0,0,1];
w10(:,:,9)= [1,0,0,0,0,0,0,0,0,0;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,0.75,0.25;
    0,0,0,0,0,0,0,0,0,1];
w10(:,:,10)= [0.75,0,0,0,0,0,0,0,0,0.25;
    0,1,0,0,0,0,0,0,0,0;
    0,0,1,0,0,0,0,0,0,0;
    0,0,0,1,0,0,0,0,0,0; 
    0,0,0,0,1,0,0,0,0,0;
    0,0,0,0,0,1,0,0,0,0;
    0,0,0,0,0,0,1,0,0,0;
    0,0,0,0,0,0,0,1,0,0;
    0,0,0,0,0,0,0,0,1,0;
    0,0,0,0,0,0,0,0,0,1];
for i=1:time_horizon_final
   w10(:,:,i+10)=w10(:,:,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



w = zeros(N,N,time_horizon_final+10,N);
w(:,:,:,1) = w1;
w(:,:,:,2) = w2;
w(:,:,:,3) = w3;
w(:,:,:,4) = w4;
w(:,:,:,5) = w5;
w(:,:,:,6) = w6;
w(:,:,:,7) = w7;
w(:,:,:,8) = w8;
w(:,:,:,9) = w9;
w(:,:,:,10) = w10;