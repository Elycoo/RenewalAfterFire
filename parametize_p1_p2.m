%% determine if sum final result of field boting of not
%% different resolution for p1 p2

p1_vec = 1:0.1:2;
len1 = length(p1_vec);
len2 = 10;
p2_vec = zeros(len1,len2);

for ii = 1:len1
    p2_vec(ii,:) = linspace(max(1,p1_vec(ii)-0.03),min(2,p1_vec(ii)+0.0),len2);
end


kBT_vec = 1.3; % kBT higher -> more changes in field

len3 = length(kBT_vec);

fprintf('Estimation time to run (ETA) is %.3g minutes.\n',len1*len2*len3*20/60)
timer_all = tic();
result = zeros(len1,len2,len3);
boring = 0;

%%

runFromOtherScript = true;
for cc = 1:len3
for bb = 1:len2
for aa = 1:len1
    
    p1 = p1_vec(aa);
    p2 = p2_vec(aa,bb);
    kBT = kBT_vec(cc);
    RenewalAfterFire;
    num = [numOfA(end) numOfB(end) numOfE(end)]/Size^2;
    if any(num > 0.9)
        % if one type is govern all other types
        result(aa,bb,cc) = boring;
        continue;
    elseif mean(diff(numOfA(end-200:end))) < -1 || mean(diff(numOfB(end-200:end))) < -1
        % if one type is in the way to govern all
        result(aa,bb,cc) = boring;
        continue;
    end
    
    analyze_blub;
    small_groups = [values_of_histA(1,end) values_of_histA(1,end)]/Size^2;
    if all(small_groups > 1000/Size^2)
        % if most of the area in the last frame is in small bulb and there
        % is no big cluster so it's random like
        result(aa,bb,cc) = boring;
        continue;
    elseif sum(mean_A(20:end)) < 0.1 && sum(mean_B(20:end)) < 0.1 % it's allowed to be about 0.5
        % if there is no big groups
        result(aa,bb,cc) = boring;
        continue;
    end
    
    % if there is no somthing of above so maybe it's not boring
    result(aa,bb,cc) = ~boring;
        
end
end
end
toc(timer_all)/60
%%
% [aa,bb,cc] = ind2sub(size(result),find(result==1));
% ii = 13;
% sum(sum(sum(result)))
% p1 = p1_vec(aa(ii));
% p2 = p2_vec(bb(ii));
% kBT = kBT_vec(cc(ii));
% 


for ii = 1:1:size(result,3)
figure(1);
imagesc(result(:,:,ii))
title(['kBT is ' num2str(kBT_vec(ii))])
colorbar;
xlabel('g - make growing')  
ylabel('r - interupt growing')
set(gca,'YDir','normal')
set(gca,'FontSize',20)
set(gcf,'Position',[334   231   800   700])
% colorbar('Ticks',[0,1],...
%          'TickLabels',{'No Islands','With Islands'})
pause(1)
end

% makeMyGif(result, 'map of boringness3.gif',1);

%%
save('results_map_of_boringness4_parametize_p1_p2_2',...
    'result','p1_vec', 'p2_vec', 'kBT_vec',...
    'Size', 'N','A_shape' ,'B_shape','E_shape',...
    'depth_of_shape','depth_of_shape_for_E','iter_to_be_stat')