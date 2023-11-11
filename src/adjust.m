function [theta, w] = adjust(T_j, gamma_j, val)
    %Get the eigenvectors and eigenvalues of the adjusted T_j
    %to compute the corresponding bounds according to article [1]


    j=size(T_j,1);

    %Solve the system to get phi
    vec1 = [zeros(j-1,1);gamma_j^2];
    sol = (T_j - val*eye(j))\vec1;
    phi = val+sol(end);

    %Adjust T_j
    vec2 = [zeros(j-1,1);gamma_j];
    adjusted_T_j = [T_j  vec2];
    adjusted_T_j = [adjusted_T_j ; vec2' phi];
    
    %Extract its eigenvalues and eigenvectors to compute the Gauss rule
    [V,D] = eig(adjusted_T_j);

    theta = diag(D);
    w = V(1,:)';

end