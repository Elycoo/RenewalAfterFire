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
A = rand(Size);
B = rand(Size);
E = rand(Size);

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

ae = 10;
be = 10;

ea = 1;
eb = 1;

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
Hea = @(x) -ea*x;
Heb = @(x) -ea*x;

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
% 0.7 + 0.2/N*x % this line make the probabilty after 80 test to be 0.7 of more
P = @(x) exp( log(0.7 + 0.25/N/2*x) /N );



[II,JJ] = meshgrid(1:Size);

numOfA = zeros(1,N);
distance_to_prev = zeros(1,N);
distance_to_start = zeros(1,N);
f_start = f;
count = 0;
count2 = 0;

for n = 1:N
    RandIndex = randperm(Size^2);
    f_prev = f;
    
for kk = RandIndex
    ii = II(kk);
    jj = JJ(kk);
    
    f_site = f(ii,jj);
    
    % check nn - nearest neighborhoods
    nn = [f(mod(ii-2,Size)+1,jj) f(mod(ii,Size)+1,jj) f(ii,mod(jj-2,Size)+1) f(ii,mod(jj,Size)+1)];
    a = sum(nn == 1);
    b = sum(nn == -1);
    e = sum(nn == 0);
    
    Ea = (tanh(Ha(a,b,e))+1)./2;
    Eb = (tanh(Hb(a,b,e))+1)./2;
    Ee = (tanh(He(a,b,e))+1)./2;
    
    
    if f_site == 1 % f_site == A
        Pa = P(n);
        Pb = (1-Pa)*Eb./(Eb+Ee);
        Pe = (1-Pa)*Ee./(Eb+Ee);
    elseif f_site == -1 % f_site == B
        Pb = P(n);
        Pa = (1-Pb)*Ea./(Ea+Ee);
        Pe = (1-Pb)*Ee./(Ea+Ee);
    elseif f_site == 0 % f_site == E
        Pe = P(n);
        Pa = (1-Pe)*Ea./(Eb+Ea);
        Pb = (1-Pe)*Eb./(Eb+Ea);
    end
%     Ea = (tanh(Ha(a,b,e))+1)./2;
%     Eb = (tanh(Hb(a,b,e))+1)./2;
%     Ee = (tanh(He(a,b,e))+1)./2;
%     
%     Ea = (atan(Ha(a,b,e))+pi./2)./pi;
%     Eb = (atan(Hb(a,b,e))+pi./2)./pi;
%     Ee = (atan(He(a,b,e))+pi./2)./pi;
%     
%     E = [Eb Ee Ea];
%     I = f_site+2;
%     I_Logic = false(1,3);
%     I_Logic(I) = true;
%     E(~I_Logic) = 0;
%     
%     self_factor = 300;
%     Etot = self_factor*E(I) + sum([Eb Ee Ea]);
%     
%     Ea = (Ea+self_factor*E(3))./Etot;
%     Eb = (Eb+self_factor*E(1))./Etot;
%     Ee = 0*(Ee+self_factor*E(2))./Etot;
% 
%     Etot = Ea + Eb + Ee;
%     Ea = Ea./Etot;
%     Eb = Eb./Etot;
%     Ee = Ee./Etot;

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
end

if mod(n,10) == 0
    figure;
    imagesc(f)
    colorbar
end

numOfA(n) = sum(sum(f==1));
distance_to_prev(n) = sum(sum(f~=f_prev));
distance_to_start(n) = sum(sum(f~=f_start));
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
distance_to_start(N)./Size.^2

%%
%
num = get(gcf,'Number');
for ii = num:-1:1
    figure(ii)
end

for ii = num-2:num
    figure(ii)
end
    
%}