function MakeAnimation(p, t, U, fileName, frameTime)
    h = figure;
    axis tight manual % this ensures that getframe() returns a consistent size
    for frameIdx = 1:size(U, 2)
        % Draw plot
        pdesurf(p, t, U(:,frameIdx))
        xlim([-4,2])
        ylim([-4,2])
        zlim([0,12])
        drawnow
        % Capture the plot as an image
        frame = getframe(h);
        im = frame2im(frame); 
        [imind,cm] = rgb2ind(im,256); 
        % Write to the GIF File
        if frameIdx == 1 
            imwrite(imind,cm,[fileName, '.gif'],'gif', 'Loopcount',inf,'DelayTime',frameTime); 
        else 
            imwrite(imind,cm,[fileName, '.gif'],'gif','WriteMode','append','DelayTime',frameTime); 
        end
    end
end