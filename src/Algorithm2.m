function [Up_list,Lp_list,I_list,Itermin,Itermax]=Algorithm2(fun, A, maxit, tol, rep, p)
%Compute the trace estimate of the matrix function fun(A) given a
%tolerance tol, a probabibilty p for the CI, a number of repetition rep and a
%maximum number of iterations maxit for Algorithm 1

%Return the list of the confidence bound as well as the List of the
%approximation as the iterations go on

    I_list=[];
    L_list=[];
    U_list=[];
    Lp_list=[];
    Up_list=[];
    n=size(A,1);

    for j=1:rep
        rng(j);
        %Compute the rademacher random vector
        z=2*(rand(n,1)<0.5)-1;

        [U, L, Iter]=Algorithm1(fun, A, z, maxit, tol);
        U_list=[U_list U];
        L_list=[L_list L];

        if j==1
            Umax=U;
            Lmin=L;
            Itermin = Iter;
            Itermax = Iter;
        else
            Umax=max(U, Umax);
            Lmin=min(L, Lmin);
            Itermin = min(Iter,Itermin);
            Itermax = max(Iter,Itermax);
        end

        %Compute the estimate and the confidence interval
        I_list=[I_list mean([U_list L_list])];
        eta=sqrt(-0.5*j*(Umax-Lmin)^2*log((1-p)/2));
        Lp_list=[Lp_list mean(L_list)-eta/j];
        Up_list=[Up_list mean(U_list)+eta/j];
    end
    I_p=I_list(end);
    L_p=Lp_list(end);
    U_p=Up_list(end);
end