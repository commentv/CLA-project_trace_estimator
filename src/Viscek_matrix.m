clear;

%%% Parameters for the experiment %%%
tol = 1e-4; %Tolerance for Algorithm 1

%%%Initial conditions to form the VFH matrix %%%
H = [-4 1 1 1 1;
    1 -2 0 0 0;
    1 0 -2 0 0;
    1 0 0 -2 0;
    1 0 0 0 -2];
H = sparse(H);
%The p_j_0 are 1 for the relation in the article to hold
p_1_k_2 = 1;
p_2_k_2 = 1;
p_3_k_2 = 1;
p_4_k_2 = 1;

for i = 1:3
    [H,p_1_k_1,p_2_k_1,p_3_k_1,p_4_k_1] = formVFH(H,p_1_k_2,p_2_k_2,p_3_k_2,p_4_k_2);

    p_1_k_2 = p_1_k_1;
    p_2_k_2 = p_2_k_1;
    p_3_k_2 = p_3_k_1;
    p_4_k_2 = p_4_k_1;

end

n=length(H);
cond = condest(H);
I = eye(n);
i_vect = [1,100,301,625];%Choice of i in the article

fprintf('---VFH fractal Hamiltonian, n =%d, cond(A) = %d---\n',n,cond)
fprintf('--- A(i,i)^-1, Gauss-Radau---\n')
fprintf(' i , iter , Lower bound L_i , Upper bound U_i , U_i - L_i  \n')

for k = 1:length(i_vect)
    i = i_vect(k);
    u = I(:,i);

    [U,L,iter] = Algorithm1(@(x) 1./x,-H,u,n,tol);
    % We took -H because H as defined in the article is negative definite,
    % To match the results, we had to take -H
    fprintf(' %d   %d     %e       %e      %e  \n',i,iter,L,U,U-L)
end


[H,p_1_k_1,p_2_k_1,p_3_k_1,p_4_k_1] = formVFH(H,p_1_k_2,p_2_k_2,p_3_k_2,p_4_k_2);

p_1_k_2 = p_1_k_1;
p_2_k_2 = p_2_k_1;
p_3_k_2 = p_3_k_1;
p_4_k_2 = p_4_k_1;

n=length(H);
cond = condest(H);
I = eye(n);
i_vect = [1,100,2000,3125];%Choice of i in the article

fprintf('---VFH fractal Hamiltonian, n =%d, cond(A) = %d---\n',n,cond)
fprintf('--- A(i,i)^-1, Gauss-Radau---\n')
fprintf(' i , iter , Lower bound L_i , Upper bound U_i , U_i - L_i  \n')

for k = 1:length(i_vect)
    i = i_vect(k);
    u = I(:,i);

    [U,L,iter] = Algorithm1(@(x) 1./x,-H,u,n,tol);
    
    fprintf(' %d   %d     %e       %e      %e  \n',i,iter,L,U,U-L)
end