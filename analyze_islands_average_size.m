function m = analyze_islands_average_size(f)
%% determine the average area of patches

A_type = f == 1; % Find A places.
labeledA = bwlabel(A_type,4);
A_area = (table2array(regionprops('table',labeledA, 'Area')));

B_type = f == -1; % Find B places.
labeledB = bwlabel(B_type,4);
B_area = (table2array(regionprops('table',labeledB, 'Area')));

m = mean([A_area;B_area]);