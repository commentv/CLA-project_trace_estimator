clear;

%%% Parameters for the experiment %%%
n = 900; %Size of the matric
tol = 1e-4; %Tolerance for algorithm 1

%%% Form the matrix and its inverse %%%
A = gallery('poisson',sqrt(n));
cond = condest(A);
I = eye(n); %To get easily the canonical vectors
i_vect = [2,1,10,41,58,450,550,600,650]; %Choice of i in the article
j_vect = [1,900,90,42,59,449,750,602,750]; %Choice of j in the article

fprintf('---Poison matrix, n =%d, cond(A) = %d---\n',n,cond)
fprintf(' i , j , iter , Lower bound L_i , Upper bound U_i , U_i - L_i  \n','Interpreter','latex')

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
    
    fprintf(' %d   %d  %d,%d    %e     %e     %e  \n',i,j,iter_y,iter_z,L,U,U-L)
end
