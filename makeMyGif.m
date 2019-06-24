function makeMyGif(all_images,filename,DelayTime,xData,yData,labelX,labelY,labelTitles)
%% Create GIF from 3D matrix of images
if DelayTime <= 0
    DelayTime = 0.1;
end
h = figure(1);
axis tight manual % this ensures that getframe() returns a consistent size
for ii = 1:size(all_images,3)
    % Draw plot for all images
    if isempty(xData)
        imagesc(all_images(:,:,ii))
%         map = [[120/255   222/255   0];
%                 [0.9769    0.9839    0.0805];
%                 [0.2422    0.1504    0.6603]];
%         title(yData(ii))
%         axis square
%         set(gca,'YDir','normal')
%         set(gcf,'Position',[575   402   665   576])
%         colormap(map)
%         colorbar('Ticks',[-2/3,0,2/3],...
%                 'TickLabels',{'B','E','A'})
%         set(gca,'FontSize',20)
%         xticks([])
%         yticks([])

    else
        imagesc(xData,yData,all_images(:,:,ii))
    end
%     title(ii)
%     colorbar;
    if iscell(labelTitles)
        title(labelTitles{ii})
    end
    if ~strcmp(labelX,'')
        xlabel(labelX)
    end
    if ~strcmp(labelY,'')
        ylabel(labelY)
    end
    set(gca,'YDir','normal')
    set(gca,'FontSize',20)
%     set(gcf,'Position',[334   231   800   700])
%     colorbar('Ticks',[0,1],...
%              'TickLabels',{'Boring','Not boring'})
%     text(1.7,1.2,'With Islands','FontSize',20)
%     text(1.0,1.8,'No Islands / Extinction','FontSize',20,'Color',[1 1 1])

    
    drawnow
    % Capture the plot as an image
    frame = getframe(h);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    % Write to the GIF File
    if ii == 1
        imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',1);
%     elseif ii == size(all_images,3)
%         imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',2);
    else
        imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',DelayTime);
    end
end


end
