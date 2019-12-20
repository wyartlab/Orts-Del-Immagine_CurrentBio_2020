function M2avi_color_ROIs(M,avi, rois);

mov = repmat(struct('cdata',[],'colormap',[]),[1 size(M,3)]);

map = jet(256);

figure(gcf);

set(gcf,'doublebuffer','on');
set(gcf,'Units','Normalized','Position', [0.25, 0.25, 0.5, 0.5]) 
writerObj=VideoWriter(avi);
open(writerObj);

couleur = [0 0 1; 0 1 0; 1 0 0 ; 1 1 0; 1 0 1; 0 1 1; 0 0 0.5; 0 0.5 0; 0.5 0 0 ; 0.5 0.5 0; 0.5 0 0.5; 0 0.5 0.5; 0 0 0.82; 0 0.82 0; 0.82 0 0 ; 0.82 0.82 0; 0.82 0 0.82; 0 0.82 0.82; 0 0 0.65; 0 0.65 0; 0.65 0 0 ; 0.65 0.65 0; 0.65 0 0.65; 0 0.65 0.65];
COLOR=[0 0 1; 0 1 0; 1 0 0];

for frame=1:size(M,3)
    
    imshow(M(:,:,frame),map);
        for i=1:size(rois,2)
            %patch('X',rois{i}(:,1),'Y',rois{i}(:,2),'EdgeColor', 'w', 'FaceColor', 'none', 'LineWidth', 1);
            patch('X',rois{i}(:,1),'Y',rois{i}(:,2),'EdgeColor', couleur(1+mod(i-1,24),:), 'FaceColor', 'none', 'LineWidth', 1);
            %text(rois{i}(1,1),rois{i}(1,2)+15,num2str(i),'Color','w', 'LineWidth', 1);
            text(rois{i}(1,1),rois{i}(1,2)+15,num2str(i),'Color',couleur(1+mod(i-1,24),:), 'LineWidth', 1);
        end;
    disp(frame)
    frame2=getframe;
    writeVideo(writerObj, frame2);
    
end

close(writerObj)

