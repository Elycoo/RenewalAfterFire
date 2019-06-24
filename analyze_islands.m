%% analyze all_f
N = size(all_f,3)-1;
edges = [1:0.5:20];

values_of_histA = zeros(numel(edges)-1,N);
values_of_histB = zeros(numel(edges)-1,N);
for ii = 1:N
    f = all_f(:,:,ii+1);
    
    A_type = f == 1; % Find A places.
    labeledA = bwlabel(A_type,4);
    A_area = sqrt(table2array(regionprops('table',labeledA, 'Area')));
    h = histcounts(A_area,edges);
    values_of_histA(:,ii) = h;
    
    B_type = f == -1; % Find B places.
    labeledB = bwlabel(B_type,4);
    B_area = sqrt(table2array(regionprops('table',labeledB, 'Area')));
    h = histcounts(B_area,edges);
    values_of_histB(:,ii) = h;
    
%     'not counted:'
%     sum(A_area>20)
%     sum(B_area>20)
end
mean_A = mean(values_of_histA,2);
mean_B = mean(values_of_histB,2);
% 
% bar(edges(1:end-1),sqrt(mean_B))
% bar(edges(1:end-1),sqrt(mean_A))