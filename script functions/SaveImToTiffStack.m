function SaveImToTiffStack(IM,tiffdir)
% IM: 1-d cell array containing all image planes

% display each plane and save as tif
h = figure;
for i_plane = 1:length(IM)
    im = IM{i_plane};
    image(im);axis equal; axis off
    drawnow;
    % save tiff
    if (i_plane == 1)
        imwrite(im, tiffdir, 'compression','none','writemode','overwrite')
    else
        imwrite(im, tiffdir, 'compression','none','writemode','append')
    end
    %     pause(0.2)
end
close(h)