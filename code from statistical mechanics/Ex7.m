%% Exercize 7,8 in statistical mechanics.
% Metropolis algorithem using monte-karlo to calculate Ising model 2D.
% (Spins in grid 20x20)

Size = 20;
Size_sq = Size^2;

% Part A J=0, H=1
J = 0;
H = 0.25;
S = randi([0 1], Size)*2-1;
%%

T = 10;
N = 400;

T_vec = exp(0.01:01:10);
T_vec = (0.001:0.1:40)/10;

Spin_avg = zeros(size(T_vec));
Spin_T = zeros(size(T_vec));
susceptibility = zeros(size(T_vec));
susceptibility2 = zeros(size(T_vec));
E_T = zeros(size(T_vec));
C_v2 = zeros(size(T_vec));
C_v = zeros(size(T_vec));

tic

t = 0; % iterate inside T loop
for T = T_vec
t = t+1;

k = 1;
beta = 1/(k*T);

S = randi([0 1], Size)*2-1;
% r = randi([1 Size] ,Size_sq*N,2);
for n = 1:N
for ii = 1:Size
for jj = 1:Size
    
    S_site = S(ii,jj);
    E0 = -H * S_site;
    S_site = - S_site;
    E1 = -H * S_site;

    DeltaE = E1 - E0;
    if DeltaE < 0
        S(ii,jj) = S_site;
    else
        if rand(1) < exp(-beta*DeltaE)
            S(ii,jj) = S_site;
        end
    end
end
end
Spin_avg(n) = sum(sum(S))/Size_sq;
end



% plot(Spin_avg)
Spin_T(t) = mean(Spin_avg(200:end));
susceptibility(t) = beta*std(Spin_avg(200:end)).^2;
C_v(t) = k*(beta.^2)*std(H*Spin_avg(200:end)).^2;


end
toc

%% plotting!
close all
hold on;
E_T = -H*Spin_T;

plot(T_vec,Spin_T,'r*');
plot(T_vec,E_T,'m*');
plot(T_vec,C_v*Size_sq,'b*');
plot(T_vec,susceptibility*Size_sq,'g*');


% analytic results

T = T_vec;
k = 1;
Bj = 2/2;
beta = 1./(k*T_vec);
x = beta*H*1.5;
Z = sinh((Bj+1/2)*x)./sinh(x/2);
M = (Bj+1/2)*coth((Bj+1/2)*x)-0.5*coth(x/2);

dH = 0.0001;
H_ = H + dH;
x_ = beta*H_*1.5;
Z_ = sinh((Bj+1/2)*x_)./sinh(x_/2);
M_ = (Bj+1/2)*coth((Bj+1/2)*x_)-0.5*coth(x_/2);


plot(T_vec,M,'k','LineWidth',2)
plot(T_vec,-H*M,'k','LineWidth',2)
plot(T_vec(1:end-1),diff(-H*M)./diff(T_vec),'k','LineWidth',2)
plot(T_vec,(M_-M)./dH,'k','LineWidth',2)



legend('<Spin>','E(T)','C_v','\chi','Analytic')

title('ISING H=0.25, J=0')
set(gca,'FontSize',30)