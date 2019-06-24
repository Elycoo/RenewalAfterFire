%% The average Size of patches

r_vec = 1:0.1:2; % parameter that interupt the growing.
g_vec = 1:0.1:2; % parameter that make growing.
kBT_vec = 1.4:0.5:1.8; % kBT higher -> more changes in field

r_vec = 1.7; % parameter that interupt the growing.
g_vec = 1.3; % parameter that make growing.
kBT_vec = 0.6:0.1:1.8; % kBT higher -> more changes in field

len1 = length(r_vec);
len2 = length(g_vec);
len3 = length(kBT_vec);

fprintf('Estimation time to run (ETA) is %.3g minuest.\n',len1*len2*len3*20/60)
timer_all = tic();
result = zeros(len1,len2,len3);
boring = 0;
last_images = zeros(100,100,len1,len2,len3);
%%

runFromOtherScript = true;
for cc = 1:len3
for bb = 1:len2
for aa = 1:len1
    
    p1 = r_vec(aa);
    p2 = g_vec(bb);
    kBT = kBT_vec(cc);
    
    RenewalAfterFire;
    last_images(:,:,aa,bb,cc) = f;
    result(aa,bb,cc) = analyze_islands_average_size(f);
        
end
end
end
toc(timer_all)/60

%%
save('results_average_size_of_islands_6',...
    'result','r_vec', 'g_vec', 'kBT_vec',...
    'Size', 'N','A_shape' ,'B_shape','E_shape',...
    'depth_of_shape','depth_of_shape_for_E','iter_to_be_stat','last_images')
%%
cell_titles = cell(size(result,3),1);
for ii = 1:size(result,3)

imagesc(r_vec,g_vec,result(:,:))
cell_titles{ii} = ['kBT is ' num2str(kBT_vec(ii))];
title(cell_titles{ii})
colorbar;
xlabel('g - make growing')
ylabel('r - interupt growing')
set(gca,'YDir','normal')
set(gca,'FontSize',20)
set(gcf,'Position',[334   231   800   700])
text(1.7,1.2,'With Islands','FontSize',20)
text(1.1,1.8,'No Islands','FontSize',20,'Color',[1 1 1])
% colorbar('Ticks',[0,1],...
%          'TickLabels',{'Boring','Not boring'})

end
plot(kBT_vec,result(:))


%% tempertue dependence
yyaxis left
plot(kBT_vec,mean(A(1:3,:),1),'.-','MarkerSize',24)
ylabel('Size of islands')
yyaxis right
plot(kBT_vec,mean(A(4:5,:),1),'.-','MarkerSize',24)
ylabel('Size of islands')


set(gca,'FontSize',20)
% set(gcf,'Position',[334   231   800   700])
xlabel('k_BT')
legend('Islands zone','Non island zone')