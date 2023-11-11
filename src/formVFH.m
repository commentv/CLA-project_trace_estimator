function [H,p_1_k_1,p_2_k_1,p_3_k_1,p_4_k_1] = formVFH(H_init,p_1_k_2,p_2_k_2,p_3_k_2,p_4_k_2)
%Perform one step to form the VFH matrix of the recurrence relation 
%To form the H_k

    power = log(length(H_init))/log(5);
    k=length(H_init);
    p_1_1 = 3;
    p_2_1 = 2;
    p_3_1 = 5;
    p_4_1 = 4;
    
    %%% Form the p_j_k_1's %%%
    p_1_k_1 = round(5^(power-1)*(p_1_1-1)+p_1_k_2);
    p_2_k_1 = round(5^(power-1)*(p_2_1-1)+p_2_k_2);
    p_3_k_1 = round(5^(power-1)*(p_3_1-1)+p_3_k_2);
    p_4_k_1 = round(5^(power-1)*(p_4_1-1)+p_4_k_2);
    
    H = blkdiag(H_init,H_init,H_init,H_init,H_init);
    H(k+p_1_k_1,p_2_k_1) = 1;
    H(p_2_k_1,k+p_1_k_1) = 1;

    H(2*k+p_2_k_1,p_1_k_1) = 1;
    H(p_1_k_1,2*k+p_2_k_1) = 1;

    H(3*k+p_3_k_1,p_4_k_1) = 1;
    H(p_4_k_1,3*k+p_3_k_1) = 1;

    H(4*k+p_4_k_1,p_3_k_1) = 1;
    H(p_3_k_1,4*k+p_4_k_1) = 1;
end