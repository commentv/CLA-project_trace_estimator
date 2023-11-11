clear;

%%% Parameters for the experiment %%%
n = 900; %Size of the matrix%
k = 30; %Size of the blocked matrices
tol = 1e-4;%Tolerance for Algorithm 1
nu = 0.2;%Parameter for Linear Heat Flow matrix

%%% Form the matrix and its inverse %%%
A = diag(-nu*ones(n-k,1),-k)+diag(-nu*ones(n-1,1),-1)+diag((1+4*nu)*ones(n,1))...
    +diag(-nu*ones(n-1,1),1)+diag(-nu*ones(n-k,1),k);
A=sparse(A);
inv_A = inv(A);
cond = condest(A);
I = eye(n); %To get easily the canonical vectors
i_vect = [1,22,32]; %Choice of i in the article

fprintf('---Linear heat flow matrix, n =%d, cond(A) = %d---\n',n,cond)
fprintf('--- A(i,i)^-1, Gauss-Radau---\n')
fprintf(' i , Exact value, iter , Lower bound L_i , Upper bound U_i , U_i - L_i  \n')

for k = 1:length(i_vect)
    i = i_vect(k);
    u = I(:,i); %Get the corresponding canonical vector

    [U,L,iter] = Algorithm1(@(x) 1./x,A,u,n,tol);
    exact = inv_A(i,i);
    
    fprintf(' %d   %e   %d     %e       %e      %e  \n',i,full(exact),iter,L,U,U-L)
end

i_vect = [2,20,200,200,899]; %Choice of i in the article
j_vect = [1,21,181,700,895]; %Choice of j in the article

fprintf('---Linear heat flow matrix, n =%d, cond(A) = %d---\n',n,cond)
fprintf('--- A(i,j)^-1, Gauss-Radau---\n')
fprintf(' i , j , Exact value, iter , Lower bound L_i , Upper bound U_i , U_i - L_i  \n','Interpreter','latex')

for k = 1:length(i_vect)
    i = i_vect(k);
    j = j_vect(k);
    u = I(:,i);
    v = I(:,j); %Get the corresponding canonical vectors
    y = u+v;
    z = u-v;

    [U_y,L_y,iter_y] = Algorithm1(@(x) 1./x,A,y,n,tol);
    [U_z,L_z,iter_z] = Algorithm1(@(x) 1./x,A,z,n,tol);

    U = (U_y-L_z)/4;
    L = (L_y-U_z)/4; %Derivation from the article
    exact = inv_A(i,j);
    
    fprintf(' %d   %d   %e  %d    %e     %e     %e  \n',i,j,full(exact),iter_y,L,U,U-L)
end
