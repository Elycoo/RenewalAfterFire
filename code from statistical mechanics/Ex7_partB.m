%% Exercize 7,8 in statistical mechanics.
% Metropolis algorithem using monte-karlo to calculate Ising model 2D.
% (Spins in grid 20x20)

Size = 20;
Size_sq = Size^2;

% Part B J=1, H=0
J = 0;
H = 0.25;
N = 400;
T_vec = (0.001:0.1:40)/10;

Spin_avg = zeros(size(T_vec));
Spin_T = zeros(size(T_vec));
susceptibility = zeros(size(T_vec));
susceptibility2 = zeros(size(T_vec));
E_T = zeros(size(T_vec));
C_v2 = zeros(size(T_vec));
C_v = zeros(size(T_vec));

%% 
tic
t = 0; % iterate inside T loop
for T = T_vec
t = t+1;
E_T(t) = 0;

k = 1;
beta = 1/(k*T);

S = randi([-1 1], Size);
S(S==0) = 1;
for n = (1:N)
for ii = 1:Size
for jj = 1:Size
    
    S_site = S(ii,jj);
    E0 = -J*(  S(mod(ii,Size)+1,jj)*S_site +   S(mod(ii-2,Size)+1,jj)*S_site ...
             + S(ii,mod(jj,Size)+1)*S_site +   S(ii,mod(jj-2,Size)+1)*S_site ) ...
         -H * S_site;
    S_site = - S_site;
    E1 = -J*(  S(mod(ii,Size)+1,jj)*S_site +   S(mod(ii-2,Size)+1,jj)*S_site ...
             + S(ii,mod(jj,Size)+1)*S_site +   S(ii,mod(jj-2,Size)+1)*S_site ) ...
         -H * S_site;
    

    DeltaE = E1 - E0;
    if DeltaE < 0
        S(ii,jj) = S_site;
    else
        if rand(1) < exp(-beta*DeltaE)
            S(ii,jj) = S_site;
        end
    end
    if n == N
        S_site = S(ii,jj);
        E_final = -J*(  S(mod(ii,Size)+1,jj)*S_site +   S(mod(ii-2,Size)+1,jj)*S_site ...
                      + S(ii,mod(jj,Size)+1)*S_site +   S(ii,mod(jj-2,Size)+1)*S_site ) ...
                  -H * S_site;
        E_T(t) = E_T(t) + E_final;
    end
end
end
Spin_avg(n) = sum(sum(S))/Size_sq;
end



% plot(Spin_avg)
Spin_T(t) = mean(Spin_avg(200:end));
susceptibility(t) = beta*std(Spin_avg(200:end)).^2;
C_v(t) = (k*beta.^2)*std(H*Spin_avg(200:end)).^2;


end %T end
toc

%% plotting!
close all
hold on;
plot(T_vec,Spin_T,'r*');
plot(T_vec,E_T/Size_sq,'m*');

plot(T_vec,C_v*Size_sq,'*b-');
plot(T_vec,susceptibility*Size_sq,'g*-')

% C_v2 = diff(E_T)./diff(T_vec);
% susceptibility2 = -diff(Spin_T)./diff(T_vec);
% plot(T_vec(1:end-1),C_v2,'b');
% plot(T_vec(1:end-1),susceptibility2,'*-g')

Bj = 1;
s = 1;
n = 2;      % Dimension of the problem
Tc = (2/3)*(s*(s+1))*n*Bj/k;

vec = Tc+0.01:0.001:4;
plot(vec,.004./(vec-Tc),'k','LineWidth',2)

legend('<Spin>','E(T)','C_v','\chi','Analytic - A/(T-T_c)')

title('ISING H=0.0012, J=1')
set(gca,'FontSize',30)