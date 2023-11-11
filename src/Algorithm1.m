function [U, L, iter] = Algorithm1(fun, A, x, maxit, tol)
% Perform m steps of the Lanczos process, storing the full basis U
% with reorthogonalization on Hermitian A with starting vector x.
%
% This code will malfunction if A is not square, if x is of different
% dimension than A, if x = 0, or if m is greater than the dimemion
% of A.
renorm = norm(x)^2;
x = x / norm(x);
%Compute the first and last eigenvalue to get later lower and upper bound
%The use of eigs function will be discussed in the report
b = eigs(A,1);
a = eigs(A,1,'smallestabs');

T_j=zeros(1,1);

gamma_j_1=0;
U_j_1 = inf;
L_j_1 = inf;
I_j_1 = inf;
x_j_2 = 0;

for j = 1:maxit
    %Compute the coefficient from the lanczos method and update T_j
    alpha_j = x'*A*x;
    T_j(end,end)=alpha_j;
    r_j = A*x-alpha_j*x-gamma_j_1*x_j_2;
    gamma_j = norm(r_j);

    %Get the eigenvectors and eigenvalues of the adjusted T_j and compute
    %the corresponding bounds
    [theta_a,w_a] = adjust(T_j,gamma_j,a);
    U_j = sum(fun(theta_a).*w_a.^2);

    [theta_b,w_b] = adjust(T_j,gamma_j,b);
    L_j = sum(fun(theta_b).*w_b.^2);
    I_j = (U_j+L_j)/2;

    if abs(U_j-U_j_1) < tol*abs(U_j) && abs(L_j-L_j_1) < tol*abs(L_j)
        break;
    end

    %Update the different recurrences
    T_j = [T_j [zeros(j-1,1); gamma_j]];
    T_j = [T_j; [zeros(1,j-1) gamma_j 0]];
    U_j_1 = U_j;
    L_j_1 = L_j;
    I_j_1 = I_j;
    gamma_j_1 = gamma_j;
    x_j_2=x;
    x = r_j/gamma_j;
end

U = U_j * renorm;
L = L_j * renorm;
iter = j;

end