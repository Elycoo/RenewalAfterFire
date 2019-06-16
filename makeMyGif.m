function makeMyGif(all_images,filename,DelayTime)
%% Create GIF from 3D matrix of images

if DelayTime <= 0
    DelayTime = 0.1;
end
h = figure(102);
axis tight manual % this ensures that getframe() returns a consistent size
for ii = 1:size(all_images,3)
    % Draw plot for all images
    imagesc(all_images(:,:,ii))
    title(ii)
    colorbar;
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if ii == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',DelayTime); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',DelayTime); 
      end 
end


end
