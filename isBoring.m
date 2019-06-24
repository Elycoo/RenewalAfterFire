%% determine if sum final result of field boting of not

r_vec = 1:0.1:2; % parameter that interupt the growing.
g_vec = 1:0.1:2; % parameter that make growing.
kBT_vec = 1.4:0.1:1.8; % kBT higher -> more changes in field

len1 = length(r_vec);
len2 = length(g_vec);
len3 = length(kBT_vec);

fprintf('Estimation time to run (ETA) is %.3g hours.\n',len1*len2*len3*20/60/60)
timer_all = tic();
result = zeros(len1,len2,len3);
boring = 0;

%%

runFromOtherScript = true;
for cc = 1:len3
for bb = 1:len2
for aa = 1:len1
    
    r = r_vec(aa);
    g = g_vec(bb);
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
    
    analyze_islands;
    small_groups = [values_of_histA(1,end) values_of_histA(1,end)]/Size^2;
    if all(small_groups > 1000/Size^2)
        % if most of the area in the last frame is in small islands and there
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
[aa,bb,cc] = ind2sub(size(result),find(result==1));
ii = 13;
sum(sum(sum(result)))
r = r_vec(aa(ii));
g = g_vec(bb(ii));
kBT = kBT_vec(cc(ii));

cell_titles = cell(size(result,3),1);

for ii = 1:1:size(result,3)
figure(1);
imagesc(result(:,:,ii))
cell_titles{ii} = ['kBT is ' num2str(kBT_vec(ii))];
title(cell_titles{ii})
% colorbar;
xlabel('g - make growing')
ylabel('r - interupt growing')
set(gca,'YDir','normal')
set(gca,'FontSize',20)
set(gcf,'Position',[334   231   800   700])
text(8,2,'With Islands','FontSize',20)
text(2,8,'No Islands','FontSize',20,'Color',[1 1 1])
% colorbar('Ticks',[0,1],...
%          'TickLabels',{'Boring','Not boring'})
pause(1)
end

makeMyGif(result, 'map of islands 4.gif',1,r_vec,g_vec,'g - make growing','r - interupt growing',cell_titles)

%%
save('results_map_of_boringness2',...
    'result','p1_vec', 'p2_vec', 'kBT_vec',...
    'Size', 'N','A_shape' ,'B_shape','E_shape',...
    'depth_of_shape','depth_of_shape_for_E','iter_to_be_stat')

% parametize_p1_p2