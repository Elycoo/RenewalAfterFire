%% Order of Magnitude in Physics problems
%% Project of Netta Yosef & Elyasaf Cohen
%% Model the renewal of forest after fire.
% We'll using Ising Model to calculate what will be more likely to grown in
% each area cell of the forest. 
%% Assumption:
% There are some initial value of probabilty to what is likely to grow in
% each area cell.
% For simplicty there will be 2 different types of plants. (A/B)
    
countinue_previous = false;
% countinue_previous = true;
if ~countinue_previous
%     clearvar -except countinue_previous
    countinue_previous = false;
end

toPlot = 0;

if runFromOtherScript
    countinue_previous = false;
    toPlot = false;
end

% clc
% close all

%% List of all the parameter in the simulation
Size = 100; % Size of square grid.
N = 250; % The number of itaraion until stop
if ~runFromOtherScript
    r = 1.7; % parameter that interupt the growing.
    g = 1.9; % parameter that make growing.
    kBT = 1; % kBT higher -> more changes in field
end

% % the shape of the potential is determine by these lines
A_shape = [4 2 0];
B_shape = [0 2 4];
E_shape = [4 0 4];
depth_of_shape = 2;
depth_of_shape_for_E = 1;
iter_to_be_stat = 100;
% for example:
% E = (depth_of_shape + f_hist/iter_to_be_stat)*A_shape;
%%
% A represent the probabilty of each point in matrix to grow an A tree.
% B represent the probabilty of each point in matrix to grow an B tree.
% E represent the probabilty of each point in matrix to extinction.

if ~countinue_previous
rng('default');
% rand(Size);
% rand(Size);
A = rand(Size,'single');
B = rand(Size,'single');
E = 0.3*rand(Size,'single');

% nomalize probabilty:
S = A+B+E;
A = A./S;
B = B./S;
E = E./S;

% Concatenate the matriecs along 3rd dimension then use max.
S = cat(3,B,E,A);
[~,f] = max(S,[],3);
f = f-2;

% arbitraty initial condition
% f(:) = -1;
% f(30:70,[54:55]) = 1;
% f(48:51,48:51) = 1;

% Save initial condtion
A = f==1;
B = f==-1;
end

all_f = zeros(Size,Size,N+1,'single');
all_f(:,:,1) = f;

% Present the initial condtion
% figure;
% imagesc(f)
% colorbar
% title(1)

% Now f represnt the field -> Each point in f matrix tell what kind of
% plants do we have at that point. 1 for A, -1 for B, 0 for E.
%
% End of initializing the system
%% The Hamiltonian
% Next, we define the "energy" that mean how the plants interact each
% other. There will be positive interaction that the tree prefer more of
% his kind.
% On the other side we'll define also negetive interaction that trees are
% gives somebenefits to other planets and surpress his kinds of plants.

% First we'll try to consider just nearest neighborhoods.

% different energies
% r = 1; % parameter that interupt the growing.
% g = 2; % parameter that make growing.
Ha = @(a,b,e) (a+b)*r - a*g; 
Hb = @(a,b,e) (a+b)*r - b*g;
He = @(a,b,e) -(a+b)*r + (a+b)*g;
H = [r-g r 0; r r-g 0; g-r g-r 0];
H = H([2 3 1],:);
%% Monte-Carlo
% N is the number of the iteration for our model we choose period of 20
% years and we assumed that nothing chages in less then a season.
rng('default');
% rand(Size);
% rand(1);

[II,JJ] = meshgrid(1:Size);

numOfA = zeros(1,N,'single');
numOfB = zeros(1,N,'single');
numOfE = zeros(1,N,'single');
EnergyOfSystem = zeros(1,N,'single');
mean_field = zeros(1,N,'single');
c_v = zeros(1,N,'single');

distance_to_prev = zeros(1,N,'single');
distance_to_start = zeros(1,N,'single');
distance_to_all_f2 = zeros(1,N,'single');
f_start = f;

if ~countinue_previous
f_history = ones(size(f));
end

tic
for n = 1:N
    RandIndex = randperm(Size^2);
    f_prev = f;
    
    randArray = rand(Size,'single');
    measureE = 0;
for kk = RandIndex
    ii = II(kk);
    jj = JJ(kk);
    
    f_original = 1*A(ii,jj)+(-1)*B(ii,jj);
    f_original = f(ii,jj);
    f_site = f(ii,jj);
    
    % check nn - nearest neighborhoods
    nn = [f(mod(ii-2,Size)+1,jj) f(mod(ii,Size)+1,jj) f(ii,mod(jj-2,Size)+1) f(ii,mod(jj,Size)+1)];
    a = sum(nn == 1);
    b = sum(nn == -1);
    e = sum(nn == 0);
    
%     Ea = Ha(a,b,e);
%     Eb = Hb(a,b,e);
%     Ee = He(a,b,e);
%     E = (H+H')*[a;b;e]/2;
    E = H*[a;b;e];
    
    f_hist = f_history(ii,jj);
    if f_original == 1 % == A
        E_change = (depth_of_shape + f_hist/iter_to_be_stat)*A_shape; % energy cost for [B E A]
    elseif f_original == -1 % == B    
        E_change = (depth_of_shape + f_hist/iter_to_be_stat)*B_shape; % energy cost for [B E A]
    elseif f_original == 0 % == E
        E_change = depth_of_shape_for_E*E_shape; % energy cost for [B E A]
    end
    E = E_change + E';
    
    
%     f_site = 0;
    I = [-1 0 1];
    I(I==f_site) = [];
    deltaE = E(I+2) - E(f_site+2);
    if any(deltaE<0)
        [~,inx] = min(deltaE);
        f(ii,jj) = I(inx);
        measureE = measureE + min(deltaE);
    else
        P = exp(-deltaE/kBT);
        P = P./(sum(P)+1); % +1 for the state itself deltaE = 0; = P./Z;
        r = randArray(ii,jj);
        if r < P(1)
            f(ii,jj) = I(1);
            measureE = measureE + deltaE(1);
        elseif r < P(1)+P(2)
            f(ii,jj) = I(2);
            measureE = measureE + deltaE(1);
        else
            f(ii,jj) = f(ii,jj);
        end
    end
    
    if f_site == f(ii,jj)
        f_history(ii,jj) = f_history(ii,jj) + 1;
    else
        f_history(ii,jj) = 1;
    end

end
% trees died after some itarations
% f(f_history > 150) = 0;


all_f(:,:,n+1) = f;

% if mod(n,5) == 0
%     figure(102);
%     imagesc(f)
%     colorbar
%     title(n)
%     drawnow;
% end

% gather information about the process
numOfA(n) = sum(sum(f==1));
numOfB(n) = sum(sum(f==-1));
numOfE(n) = sum(sum(f == 0));
EnergyOfSystem(n+1) = EnergyOfSystem(n)+measureE;
mean_field(n) = mean(f(:));
c_v(n) = std(f(:));

distance_to_prev(n) = sum(sum(f~=f_prev));
distance_to_start(n) = sum(sum(f~=f_start));
distance_to_all_f2(n) = sum(sum(f~=all_f(:,:,2)));

end
timer = toc;
fprintf('Time for each iter is %.3g milli seconds\nTotal time is %.3g seconds\n',1000*timer/N,timer);
%% plot intersting graph that help to see if it converge and what the rate of chanses is
if toPlot
figure;
plot(numOfA)
title('numOfA')
figure;
bar(sqrt(histcounts(distance_to_prev)))
% histogram(distance_to_prev)
title('distance from prev f')
figure;
plot(distance_to_start)
title('distance from start f')
sumOfE = sum(numOfE)
distance_to_start(N)./Size.^2
distance_to_all_f2(N)./Size.^2
figure;
plot(EnergyOfSystem)
title('Energy Of System')
% figure;
% plot(mean_field)
% title('mean_field')
% figure;
% plot(c_v)
% title('c_v')

end

%% show all fields
%
for ii = 1:1:size(all_f,3)
    figure(102);
    imagesc(all_f(:,:,ii))
    title(ii)
    colorbar;
    caxis([-1 1])
%     pause(0.1)
    
end
%}

% makeMyGif(all_f(:,:,:),'evoving_r_1.2_g_1.8_T_1.4__3.gif',0.05,[],[],'','','')
