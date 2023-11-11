clear;

%%% Parameters for the experiment %%%
tol_vect = [1e-2,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8]; %Tolerances used
n = 900; %Size of the matrix
p = 0.95; %Probability of the CI in algo 2
rep = 50; %Number of iterations in algo 2
repet = 10; %For every tolerance, we repeat repet times to average the measures

Poisson_mat = gallery('poisson',sqrt(n));
I = eye(n); %To compute easily the canonical vectors

tr_inv_time_vect = [];
alg_1_time_vect = [];
lin_syst_time_vect = [];
alg_2_time_vect = [];

for tol = tol_vect

    tr_inv_time =0;
    alg_1_time = 0;
    lin_syst_time = 0;
    alg_2_time = 0;

    for j = 1:repet
        %%% Compute directly tr(inv(A)) %%%
        tic;
        tr_inv = trace(inv(Poisson_mat));
        t = toc;
        tr_inv_time = tr_inv_time +t;
        
        %%% Compute tr(inv(A)) using algorithm 1 %%%
        tic;
        tr_inv_alg_1 = 0;
        for i = 1:length(Poisson_mat)
            u = I(:,i);
            [U,L,~] = Algorithm1(@(x) 1./x,Poisson_mat,u,n,tol);
            tr_inv_alg_1 = tr_inv_alg_1 + (U+L)/2;
        end
        t = toc;
        alg_1_time = alg_1_time +t;
        
        %%% Compute tr(inv(A)) using solving n linear systems %%%
        tic;
        tr_inv_lin = 0;
        for i = 1:length(Poisson_mat)
            u = I(:,i);
            x = Poisson_mat\u;
            tr_inv_alg_1 = tr_inv_alg_1 + x(i);
        end
        t=toc;
        lin_syst_time = lin_syst_time +t;
        
        %%% Compute tr(inv(A)) using Algorithm 2 %%%
        tic;
        [~,~,tr_inv_alg_2]=Algorithm2(@(x) 1./x, Poisson_mat, length(Poisson_mat), tol, rep, p);
        t=toc;
        alg_2_time = alg_2_time + t;
   
    end

    tr_inv_time_vect = [tr_inv_time_vect tr_inv_time/repet];
    alg_1_time_vect = [alg_1_time_vect alg_1_time/repet];
    lin_syst_time_vect = [lin_syst_time_vect lin_syst_time/repet];
    alg_2_time_vect = [alg_2_time_vect alg_2_time/repet];

end
figure;
ax_1 = subplot(1,1,1,'XScale', 'log', 'YScale', 'log');
title(ax_1,'Comparison of time execution')
ylabel(ax_1,'Time(s)')
xlabel(ax_1,'Tolerance');
set(ax_1,'Xdir','reverse');
hold(ax_1,'on')
loglog(ax_1,tol_vect,tr_inv_time_vect,'b');
loglog(ax_1,tol_vect,alg_1_time_vect,'r');
loglog(ax_1,tol_vect,alg_2_time_vect,'k');
loglog(ax_1,tol_vect,lin_syst_time_vect,'m');
legend(ax_1,'tr(inv(A))','Algorithm 1','Algorithm 2','n linear systems','Location','best');
hold(ax_1,'off')