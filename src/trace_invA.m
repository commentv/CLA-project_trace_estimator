clear; 

%%% Parameters for the experiment %%%
tol = 1e-4;%Tolerance used in Algorithm 2
rep = 50;%Number of repetitions in Algorithm 2
p = 0.95;%Probability for the CI in Algorithm 2
rng(1);

%%% First, we compute the graph %%%
%%% Poisson matrix %%%
Poisson_mat = gallery('poisson',30);
cond_Poisson=condest(Poisson_mat);
[U_Poisson_list,L_Poisson_list,I_Poisson_list,IterPoissonMin,IterPoissonMax]...
    =Algorithm2(@(x) 1./x, Poisson_mat, length(Poisson_mat), tol, rep, p);
inv_Poisson_mat = inv(Poisson_mat);
exact_Poisson_vect = trace(inv_Poisson_mat) * ones(1,length(U_Poisson_list));
rel_Poisson = 100*abs((exact_Poisson_vect(end)-I_Poisson_list(end))/exact_Poisson_vect(end));
figure;
ax_1 = subplot(1,1,1,'XScale', 'linear', 'YScale', 'linear');
title(ax_1,'Poisson matrix (n=900)')
ylabel(ax_1,'Trace')
xlabel(ax_1,'Number of samples');
hold(ax_1,'on')
plot(ax_1,1:length(U_Poisson_list),U_Poisson_list,'--b');
plot(ax_1,1:length(L_Poisson_list),L_Poisson_list,'--r');
plot(ax_1,1:length(I_Poisson_list),I_Poisson_list,'-xm');
plot(ax_1,1:length(exact_Poisson_vect),exact_Poisson_vect,'-k');
legend(ax_1,'Upper CI','Lower CI','Estimated value','Exact value','Location','best');
hold(ax_1,'off')

%%% Wathen matrix %%%
Wathen_mat = gallery('wathen',12,12);
cond_Wathen=condest(Wathen_mat);
[U_Wathen_list,L_Wathen_list,I_Wathen_list,IterWathenMin,IterWathenMax]...
    =Algorithm2(@(x) 1./x, Wathen_mat, length(Wathen_mat), tol, rep, p);
inv_Wathen_mat = inv(Wathen_mat);
exact_Wathen_vect = trace(inv_Wathen_mat) * ones(1,length(U_Wathen_list));
rel_Wathen = 100*abs((exact_Wathen_vect(end)-I_Wathen_list(end))/exact_Wathen_vect(end));
figure;
ax_2 = subplot(1,1,1,'XScale', 'linear', 'YScale', 'linear');
title(ax_2,'Wathen matrix (n=481)')
ylabel(ax_2,'Trace')
xlabel(ax_2,'Number of samples');
hold(ax_2,'on')
plot(ax_2,1:length(U_Wathen_list),U_Wathen_list,'--b');
plot(ax_2,1:length(L_Wathen_list),L_Wathen_list,'--r');
plot(ax_2,1:length(I_Wathen_list),I_Wathen_list,'-xm');
plot(ax_2,1:length(exact_Wathen_vect),exact_Wathen_vect,'-k');
legend(ax_2,'Upper CI','Lower CI','Estimated value','Exact value','Location','best');
hold(ax_2,'off')

%%% Lehmer matrix %%%
Lehmer_mat = gallery('lehmer',200);
cond_Lehmer=condest(Lehmer_mat);
[U_Lehmer_list,L_Lehmer_list,I_Lehmer_list,IterLehmerMin,IterLehmerMax]...
    =Algorithm2(@(x) 1./x, Lehmer_mat, length(Lehmer_mat), tol, rep, p);
inv_Lehmer_mat = inv(Lehmer_mat);
exact_Lehmer_vect = trace(inv_Lehmer_mat) * ones(1,length(U_Lehmer_list));
rel_Lehmer = 100*abs((exact_Lehmer_vect(end)-I_Lehmer_list(end))/exact_Lehmer_vect(end));
figure;
ax_3 = subplot(1,1,1,'XScale', 'linear', 'YScale', 'linear');
title(ax_3,'Lehmer matrix (n=200)')
ylabel(ax_3,'Trace')
xlabel(ax_3,'Number of samples');
hold(ax_3,'on')
plot(ax_3,1:length(U_Lehmer_list),U_Lehmer_list,'--b');
plot(ax_3,1:length(L_Lehmer_list),L_Lehmer_list,'--r');
plot(ax_3,1:length(I_Lehmer_list),I_Lehmer_list,'-xm');
plot(ax_3,1:length(exact_Lehmer_vect),exact_Lehmer_vect,'-k');
legend(ax_3,'Upper CI','Lower CI','Estimated value','Exact value','Location','best');
hold(ax_3,'off')

%%% VFH matrix %%%
VFH_mat = [-4 1 1 1 1;
    1 -2 0 0 0;
    1 0 -2 0 0;
    1 0 0 -2 0;
    1 0 0 0 -2];
VFH_mat = sparse(VFH_mat);
p_1_k_2 = 1;
p_2_k_2 = 1;
p_3_k_2 = 1;
p_4_k_2 = 1;

for i = 1:3
    [VFH_mat,p_1_k_1,p_2_k_1,p_3_k_1,p_4_k_1] = formVFH(VFH_mat,p_1_k_2,p_2_k_2,p_3_k_2,p_4_k_2);

    p_1_k_2 = p_1_k_1;
    p_2_k_2 = p_2_k_1;
    p_3_k_2 = p_3_k_1;
    p_4_k_2 = p_4_k_1;

end
VFH_mat = -VFH_mat;
cond_VFH=condest(VFH_mat);
[U_VFH_list,L_VFH_list,I_VFH_list,IterVFHMin,IterVFHMax]...
    =Algorithm2(@(x) 1./x, VFH_mat, length(VFH_mat), tol, rep, p);
inv_VFH_mat = inv(VFH_mat);
exact_VFH_vect = trace(inv_VFH_mat) * ones(1,length(U_VFH_list));
rel_VFH = 100*abs((exact_VFH_vect(end)-I_VFH_list(end))/exact_VFH_vect(end));
figure;
ax_4 = subplot(1,1,1,'XScale', 'linear', 'YScale', 'linear');
title(ax_4,'VFH matrix (n=625)')
ylabel(ax_4,'Trace')
xlabel(ax_4,'Number of samples');
hold(ax_4,'on')
plot(ax_4,1:length(U_VFH_list),U_VFH_list,'--b');
plot(ax_4,1:length(L_VFH_list),L_VFH_list,'--r');
plot(ax_4,1:length(I_VFH_list),I_VFH_list,'-xm');
plot(ax_4,1:length(exact_VFH_vect),exact_VFH_vect,'-k');
legend(ax_4,'Upper CI','Lower CI','Estimated value','Exact value','Location','best');
hold(ax_4,'off')


%%% We then compute the table %%%
fprintf('--- A(i,i)^-1, Gauss-Radau---\n')
fprintf(' Matrix  , n   , Exact        , Iter , Estimated    , Rel. err, Confidence bounds  \n')
fprintf(' Poisson , %d , %e , %d-%d , %e , %.2f%%, (%e,%e)  \n',length(Poisson_mat),...
    exact_Poisson_vect(end),IterPoissonMin,IterPoissonMax,I_Poisson_list(end),...
    rel_Poisson,L_Poisson_list(end), U_Poisson_list(end));
fprintf(' VFH     , %d , %e , %d-%d , %e , %.2f%%, (%e,%e)  \n',length(VFH_mat),...
    exact_VFH_vect(end), IterVFHMin,IterVFHMax,I_VFH_list(end),rel_VFH,...
    L_VFH_list(end), U_Poisson_list(end));
fprintf(' Wathen  , %d , %e , %d-%d , %e , %.2f%%, (%e,%e)  \n',length(Wathen_mat),...
    exact_Wathen_vect(end),IterWathenMin,IterWathenMax,I_Wathen_list(end)...
    ,rel_Poisson,L_Wathen_list(end), U_Wathen_list(end));
fprintf(' Lehmer  , %d , %e ,  %d-%d , %e , %.2f%%, (%e,%e)  \n',length(Lehmer_mat),...
    exact_Lehmer_vect(end),IterLehmerMin,IterLehmerMax,I_Lehmer_list(end)...
    ,rel_Lehmer,L_Lehmer_list(end), U_Lehmer_list(end));