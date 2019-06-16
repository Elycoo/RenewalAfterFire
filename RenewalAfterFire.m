%% Order of Magnitude in Physics problems
%% Project of Netta Yosef & Elyasaf Cohen
%% Model the renewal of forest after fire.
% We'll using Ising Model to calculate what will be more likely to grown in
% each area cell of the forest. 
%% Assumption:
% There are some initial value of probabilty to what is likely to grow in
% each area cell.
% For simplicty there will be 2 different types of plants. (A/B)

clear
% clc
close all

%% List of all the parameter in the simulation
Size = 20; % Size of square grid.
N = 80; % The number of itaraion until stop
p1 = 1; % parameter that interupt the growing.
p2 = 2; % parameter that make growing.
beta = 0.5; % beta smaller -> more changes in field


% % the shape of the potential is determine by these lines
A_shape = [4 2 0];
B_shape = [0 2 4];
E_shape = [4 0 4];
depth_of_shape = 1.25;
depth_of_shape_for_E = 1;
iter_to_be_stat = 25;
% for example:
% E = (depth_of_shape + f_hist/iter_to_be_stat)*A_shape;

% if f_original == 1 % == A
%     E = (1.25+1*(-0.25+f_history(ii,jj))/10)*[4 2 0]; % energy cost for [B E A]
%     E_start = E;
%     E = E + [Eb Ee Ea];
% elseif f_original == -1 % == B
%     E = (1.25+1*(-0.25+f_history(ii,jj))/10)*[0 2 4]; % energy cost for [B E A]
%     E_start = E;
%     E = E + [Eb Ee Ea];
% elseif f_original == 0 % == e
%     E = 1*[4 0 4]; % energy cost for [B E A]
%     E_start = E;
%     E = E + [Eb Ee Ea];
% end



%%

% A represent the probabilty of each point in matrix to grow an A tree.
% B represent the probabilty of each point in matrix to grow an B tree.
% E represent the probabilty of each point in matrix to extinction.
% Size = 20;
rng('default');
% rand(Size);
rand(Size);
A = rand(Size);
B = rand(Size);
E = 0.3*rand(Size);

% nomalize probabilty:
S = A+B+E;
A = A./S;
B = B./S;
E = E./S;

% Concatenate the matriecs along 3rd dimension then use max.
S = cat(3,B,E,A);
[~,f] = max(S,[],3);
f = f-2;

% Save initial condtion
A = f==1;
B = f==-1;

all_f = zeros(Size,Size,N+1,'int8');
all_f(:,:,1) = f;

figure;
imagesc(f)
colorbar
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
% p1 = 1; % parameter that interupt the growing.
% p2 = 2; % parameter that make growing.
Ha = @(a,b,e) (a+b)*p1 - a*p2; 
Hb = @(a,b,e) (a+b)*p1 - b*p2;
He = @(a,b,e) -(a+b)*p1 + (a+b)*p2;
%% Monte-Carlo
% N is the number of the iteration for our model we choose period of 20
% years and we assumed that nothing chages in less then a season.
rng('default');
rand(3);
% N = 20*4;

[II,JJ] = meshgrid(1:Size);

numOfA = zeros(1,N);
numOfE = zeros(1,N);
distance_to_prev = zeros(1,N);
distance_to_start = zeros(1,N);
distance_to_all_f2 = zeros(1,N);
f_start = f;
f_history = ones(size(f));
count = 0;
count2 = 0;

for n = 1:N
    RandIndex = randperm(Size^2);
    f_prev = f;
    
    
for kk = RandIndex
    ii = II(kk);
    jj = JJ(kk);
    
    f_original = 1*A(ii,jj)+(-1)*B(ii,jj);
    f_site = f(ii,jj);
    
    % check nn - nearest neighborhoods
    nn = [f(mod(ii-2,Size)+1,jj) f(mod(ii,Size)+1,jj) f(ii,mod(jj-2,Size)+1) f(ii,mod(jj,Size)+1)];
    a = sum(nn == 1);
    b = sum(nn == -1);
    e = sum(nn == 0);
    
    Ea = Ha(a,b,e);
    Eb = Hb(a,b,e);
    Ee = He(a,b,e);
    
    f_hist = f_history(ii,jj);
    if f_original == 1 % == A
        E_start = (depth_of_shape + f_hist/iter_to_be_stat)*A_shape; % energy cost for [B E A]
    elseif f_original == -1 % == B    
        E_start = (depth_of_shape + f_hist/iter_to_be_stat)*B_shape; % energy cost for [B E A]
    elseif f_original == 0 % == e    
        E_start = depth_of_shape_for_E*E_shape; % energy cost for [B E A]
    end
    E = E_start + [Eb Ee Ea];
    
    
%     f_site = 0;
    I = [-1 0 1];
    I(I==f_site) = [];
    deltaE = E(I+2) - E(f_site+2);
    if any(deltaE<0)
        [~,inx] = min(deltaE);
        f(ii,jj) = I(inx);
    else
%         beta = 0.5; % beta smaller -> more changes in field
        P = exp(-beta*deltaE);
        P = P./(sum(P)+1); % +1 for the state itself deltaE = 0;
        r = rand(1);
        if r < P(1)
            f(ii,jj) = I(1);
        elseif r < P(2)
            f(ii,jj) = I(2);
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
all_f(:,:,n+1) = f;

if mod(n,10) == 0
    figure;
    imagesc(f)
    colorbar
    title(n)
end

numOfA(n) = sum(sum(f==1));
distance_to_prev(n) = sum(sum(f~=f_prev));
distance_to_start(n) = sum(sum(f~=f_start));
distance_to_all_f2(n) = sum(sum(f~=all_f(:,:,2)));
numOfE(n) = sum(sum(f == 0));
end
%%
figure;
plot(numOfA)
title('numOfA')
figure;
hist(distance_to_prev)
title('distance from prev f')
figure;
plot(distance_to_start)
title('distance from start f')
count
count2
numOfE = sum(numOfE)
distance_to_start(N)./Size.^2
distance_to_all_f2(N)./Size.^2

%%
num = floor(N/10)+4;
% num = get(gcf,'Number');
for ii = num:-1:1
    figure(ii)
end

for ii = num-2:num
    figure(ii)
end
    
%}

%% show all fields

for ii = 1:size(all_f,3)
%     figure(102);
%     title(ii)
%     imagesc(all_f(:,:,ii))
%     colorbar;
%     pause(0.1)
    
end
% makeMyGif(all_f, 'testAnimated3.gif',0.05);
