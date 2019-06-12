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


A = [0.2760    0.4984    0.7513    0.9593    0.8407;
     0.6797    0.9597    0.2551    0.5472    0.2543;
     0.6551    0.3404    0.5060    0.1386    0.8143;
     0.1626    0.5853    0.6991    0.1493    0.2435;
     0.1190    0.2238    0.8909    0.2575    0.9293];
% A reresent the probabilty of each point in matrix to grow an A tree.

B = [0.3500    0.3517    0.2858    0.0759    0.1299;
     0.1966    0.8308    0.7572    0.0540    0.5688;
     0.2511    0.5853    0.7537    0.5308    0.4694;
     0.6160    0.5497    0.3804    0.7792    0.0119;
     0.4733    0.9172    0.5678    0.9340    0.3371];
% B reresent the probabilty of each point in matrix to grow an B tree.
 
E = [0.1622    0.6020    0.4505    0.8258    0.1067;
    0.7943    0.2630    0.0838    0.5383    0.9619;
    0.3112    0.6541    0.2290    0.9961    0.0046;
    0.5285    0.6892    0.9133    0.0782    0.7749;
    0.1656    0.7482    0.1524    0.4427    0.8173];
% E reresent the probabilty of each point in matrix to extinction.
Size = 20;
rng('default');
rand(Size);
A = rand(Size);
B = rand(Size);
E = 0*rand(Size);

% A = 3*ones(Size);

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

all_f = zeros(Size,Size,81,'int8');
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
% n stand for the number of neighborhoods
n = 4;

% define arbitrary parameters
aa = 1;
bb = 1;
ee = 0;

ab1 = 1;
ab2 = 0;
ba1 = 1;
ba2 = 0;

ae = 1;
be = 1;

ea = -1;
eb = -1;

% A-A and B-B interaction
Haa = @(x) -aa*x.*(x-n);
Hbb = @(x) -bb*x.*(x-n);
Hee = @(x) 0;

% A-B inteaction
Hab1 = @(x) ab2*n/2;
Hab2 = @(x) ab2*n/2;
Hba1 = @(x) ba2*n/2;
Hba2 = @(x) ba2*n/2;

% A-E and A-B interaction. How A/B effected from E that around
Hae = @(x) ae*x;
Hbe = @(x) be*x;

% E-A and B-E interaction. How E effected from A/B that around
Hea = @(x) ea*x;
Heb = @(x) ea*x;

Ha = @(a,b,e) Haa(a) + Hab1(b) + Hab2(b) + Hae(e);
Hb = @(a,b,e) Hbb(b) + Hba1(a) + Hba2(a) + Hbe(e);
He = @(a,b,e) Hee(e) + Hea(a) + Heb(b);

H = @(X,a,b,e) Ha(a,b,e)*(X==1) + Hb(a,b,e)*(X==-1) + He(a,b,e)*(X==0);



%% Monte-Carlo
% N is the number of the iteration for our model we choose period of 20
% years and we assumed that nothing chages in less then a season.
rng('default');
rand(3);
N = 20*4;
% 0.7 + 0.2/N*x % this line make the probabilty after N/2 test to be 0.7 of more
P = @(x) exp( log(0.7 + 0.25/N/2*x) /(N/2) );


[II,JJ] = meshgrid(1:Size);

numOfA = zeros(1,N);
numOfE = zeros(1,N);
distance_to_prev = zeros(1,N);
distance_to_start = zeros(1,N);
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
    
    % different sigmoind function that make probabilty from x in (-Inf,Inf)
    % Option 1: not good because it decay too fast
    % Ea = (tanh(Ha(a,b,e))+1)./2;
    % Eb = (tanh(Hb(a,b,e))+1)./2;
    % Ee = (tanh(He(a,b,e))+1)./2;
    
    % Option 2: decay fine
    Ea = (atan(Ha(a,b,e))+pi./2)./pi;
    Eb = (atan(Hb(a,b,e))+pi./2)./pi;
    Ee = (atan(He(a,b,e))+pi./2)./pi;
    %     [Ea Eb Ee]
    
    if f_site == 1 % f_site == A
        if f_original == 1 || f_history(ii,jj) > 8
            Pa = P(n);
            Pb = (1-Pa)*Eb./(Eb+Ee);
            Pe = (1-Pa)*Ee./(Eb+Ee);
        else
            Pa = P(0);
            Pb = (1-Pa)*Eb./(Eb+Ee);
            Pe = (1-Pa)*Ee./(Eb+Ee);
        end
    elseif f_site == -1 % f_site == B
        if f_original == 1 || f_history(ii,jj) > 8
            Pb = P(n);
            Pa = (1-Pb)*Ea./(Ea+Ee);
            Pe = (1-Pb)*Ee./(Ea+Ee);
        else
            Pb = P(0);
            Pa = (1-Pb)*Ea./(Ea+Ee);
            Pe = (1-Pb)*Ee./(Ea+Ee);
        end
    elseif f_site == 0 % f_site == E
        if f_original == 1 || f_history(ii,jj) > 8
            Pe = 0.9*P(n);
            Pa = (1-Pe)*Ea./(Eb+Ea);
            Pb = (1-Pe)*Eb./(Eb+Ea);
        else
            Pe = 0.9*P(0);
            Pa = (1-Pe)*Ea./(Eb+Ea);
            Pb = (1-Pe)*Eb./(Eb+Ea);
        end
    end

    [~,I] = max([Pb Pe Pa]);
    I = I-2;
    
    if I ~= f_site
        f(ii,jj) = I;
        count = count + 1;
    else
        count2 = count2 + 1;
        r = rand(1);
        if r < Pa
            f(ii,jj) = 1;
        elseif r < Pa + Pe
            f(ii,jj) = 0;
        else
            f(ii,jj) = -1;
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
    figure(102);
    imagesc(all_f(:,:,ii))
    title(ii)
    colorbar;
    
    pause(0.1)
    
end
% makeMyGif(all_f, 'testAnimated2.gif',0.05);