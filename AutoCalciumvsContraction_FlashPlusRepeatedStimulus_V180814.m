function output = AutoCalciumvsContraction_FlashPlusRepeatedStimulus_V180814(framesel,ref,SelectRef);


addpath('C:')
close all

output=struc('Go', [], 'Gr', [], 'rois', [], 'dff_G', [], 'raw_G', [], 'template',[], 'contractions', []);

 %framesel= 1:total number of frames;
 %ref=image chosen for registration;
 %SelectRef=1; % If 1, selects a reference image calculated outside Matlab. If 0, calculates the reference from raw image (STD)

 
 TimeFlashContra=X % X=Time Between Flash And Contraction

 %%%%%%%%%%%%% FIRST load the movie from a multitiff file
 
%go to the right folder to perform the analysis
%set the directory for the GREEN channel (GCAMP)
disp('Choose green canal movie')
[movieG folder] = uigetfile('*.*');
G_push = strcat(folder, movieG);
Go = multitiff2M([folder movieG], framesel);

% image created to select ROIs/cells
if SelectRef==1
disp('Choose baseline STD image made on ImageJ to circle ROIs')
[image folder3] = uigetfile('*.*');
 quietG_push = multitiff2M([folder3 image],1);

else
    quietG_push = squeeze(std(single(Go),0,3));

end

fig = figure(1);
imshow(quietG_push)
set(gcf,'Units','normalized','OuterPosition',[0.1 0.1 0.8 0.8])
axis equal

pause

rois = getROIcell;


% in order to save a trace of the file used to draw ROIs and the number of each ROI 
     saveas(fig,strcat(folder,['ROI-capture-' datestr(clock,30)],'.png'));
     delete(findobj('FaceAlpha',0.49));
     delete(findobj('Color','g'));
     
%%%%%%%%%%%% SECOND REGISTER THE MOVIE TO ALIGN ALL FRAMES TO THE FIRST REFERENCE
% FRAME -- FYI A TRANSLATION in X AND Y ONLY IS USED HERE
% Click to define here the rectangle where your cells of interest are located [xmin ymin width height]
% need to check together how you load your movie and call properly the
% movie G or G_push
Answer=0;


while Answer==0
disp('now it s time to pick a SMALL rectangle in the center of the image to perform the registration')
fig = figure(1);
hold off
imagesc(quietG_push);
axis equal
rect = floor(getrect(fig));
rectindX = [rect(1):(rect(1)+rect(3))];
rectindY = [rect(2):(rect(2)+rect(4))];
template = Go(rectindY, rectindX, ref);

% bigrectX = [max(1,rect(1)-20):min(size(M,2),(rect(1)+rect(3)+20))];
% bigrectY = [max(1,rect(2)-20):min(size(M,1),(rect(2)+rect(4)+20))];

bigrectX=1:size(Go,2);
bigrectY=1:size(Go,1);

% initialize the registered movie Gr to be registered and the offset
Gr = Go;
Xoffset=[];
Yoffset=[];
TransXY=zeros(length(framesel),2);
h=waitbar(0,'Percentage of processed frames'); 
for i=1:size(Go,3)
    waitbar(i/size(Go,3),h)
    out = normxcorr2(template,Go(bigrectY,bigrectX,i));
    imshow(out);
    size(out);
    [~,ind] = max(abs(out(:)));
    [Yp,Xp] = ind2sub(size(out),ind);
    % bottom right corner matters
    Yoffset = (Yp-size(template,1))-(rect(2))+1;
    Xoffset = (Xp-size(template,2))-(rect(1))+1;
    
    TransXY(i,:)=[Yoffset,Xoffset];
    % [BW,xi,yi] = roipoly(M(:,:,1));
    se = translate(strel(1), [-Yoffset -Xoffset]);
    Gr(:,:,i) = imdilate(Go(:,:,i),se);
end
close(h)

%close all

%create avi file and check ROIs
M2avi_color_ROIs(Gr,'Gr_avi', rois);
Answer=input('is it ok?''1 if yes, 0 if no\n');
end

%Correction to avoid having -Inf after motion correction
Gr(Gr<0)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Create mask with all pixels inside ROIS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FluoAfterReg=Gr;


TotalMask=zeros(size(FluoAfterReg,1),size(FluoAfterReg,2));
Color=linspecer(size(rois,2));
MaskToPlot=TotalMask;

for i=1:size(rois,2)
    Mask{i} = roipoly(FluoAfterReg,rois{i}(:,1),rois{i}(:,2));
    TotalMask=TotalMask+Mask{i};
    MaskToPlot=MaskToPlot+double(Mask{i})*i;
    
end
TotalMask=TotalMask>0;        

%Plot ROI Images


figure(1)
hold on
imagesc(MaskToPlot)
colormap(Color)
title('InitialMask')
for i=1:size(rois,2)
    text(rois{i}(1,1),rois{i}(1,2),num2str(i),'Color','r','FontSize',16)
end
set(gca,'Ydir','reverse')


set(gcf,'Units','Normalized','Position', [0.05, 0.05, 0.9, 0.85])
saveas(gcf,strcat(folder,movieG(1:end-4),'_ROI_Pos.png'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%From the mask, calculates sum of signal inside or outside the mask
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the sum of pixels at each time inside all the ROIS

TotalActivityROIS=FluoAfterReg(repmat((TotalMask==1),1,1,size(FluoAfterReg,3)));
TotalActivityROIS=reshape(TotalActivityROIS,sum(sum(TotalMask)),size(FluoAfterReg,3));
TotalActivityROIS=sum(TotalActivityROIS)';


TotalNoise=FluoAfterReg(repmat((TotalMask==0),1,1,size(FluoAfterReg,3)));
TotalNoise=reshape(TotalNoise,sum(sum(TotalMask==0)),size(FluoAfterReg,3));
TotalNoise=sum(TotalNoise,'omitnan')';

     


%Find frames with motion artifacts and replace them by nan in Gr-->Grnan

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find Flashes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Blob=squeeze(sum(sum(FluoAfterReg>254)));
[~,locs,w,~]=findpeaks(Blob,'MinPeakHeight',5000);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Second peaks detection method : From the stimulus repetition frequency ,
%finds the max peak and deduces position of all other peaks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PeaksPos=zeros(1,length(framesel));


if (locs(end)+TimeFlashContra)>length(framesel)
    locs(end)=[];
end


NPeaks=length(locs);
 
PosFirst=locs(1)+ TimeFlashContra;
PeaksPos(locs+ TimeFlashContra)=1;
PeaksPos=PeaksPos(framesel);

figure(30)
plot(TotalActivityROIS/max(TotalActivityROIS),'linewidth',2)

hold on
plot(PeaksPos,'linewidth',2)
title('Position and widths of artifacts over all ROIS And Corrected Signal')
xlabel('Time (Frame)')
ylabel('-(Is Artifact) ')
legend('Sum All ROIS Signal','Artifact detection','Corrected Signal','Location','southeast')
set(gca,'FontSize',20,'FontWeight','bold')
set(gcf,'Units','Normalized','Position', [0.05, 0.05, 0.9, 0.85])
saveas(gcf,strcat(folder,movieG(1:end-4),'Artifact.png'));
  

%Calculates signal after contraction
PostContractionActivity=zeros(size(rois,2),NPeaks); % One baseline for the all trace
PostContractionActivity2=zeros(size(rois,2),NPeaks);% Baseline just before peak. mean on 3 points around max
PostContractionActivity3=zeros(size(rois,2),NPeaks);% Baseline just before peak. Peak = max value
Mr= Gr;
PrePeaksIdx=zeros(1,length(framesel));
for i=1:NPeaks
    PrePeaksIdx(locs(i)+ TimeFlashContra-5:locs(i)+ TimeFlashContra-3)=1;
end

raw_G=zeros(size(rois,2),length(framesel));
dff_G=zeros(size(rois,2),length(framesel));
BaselineROI=zeros(size(rois,2),1);
figure(42)
hold on
plot(PeaksPos*10-50,'linewidth',2)
Max=0;
for i=1:size(rois,2)   %For each ROI, I calculate the (spatial) average fluo signal
    bgf=Mr(repmat((Mask{i}==1),1,1,size(Mr,3)));
    bgf=reshape(bgf,sum(sum(Mask{i}==1)),size(Mr,3));
    bgf=mean(bgf)';
    
 %Calculates average Baseline on all Trace
 Baseline=mean(bgf(PrePeaksIdx==1));
 plot(Max+(bgf/Baseline -1)*100)
    
 raw_G(i,:)=bgf;
 dff_G(i,:)=((bgf/Baseline)-1)*100;
 BaselineROI(i)=Baseline;
 
 for jj=1:NPeaks
     [MaxP,PosMP]=max(bgf(3+locs(jj)+ TimeFlashContra:50+locs(jj)+ TimeFlashContra));
     Fluo=mean(bgf((PosMP+locs(jj)+ TimeFlashContra):PosMP+2+locs(jj)+ TimeFlashContra));
     Baseline2=mean(bgf(-5+locs(jj)+ TimeFlashContra:-3+locs(jj)+ TimeFlashContra));
     PostContractionActivity(i,jj)=(Fluo/Baseline-1)*100;
     PostContractionActivity2(i,jj)=(Fluo/Baseline2-1)*100;
     PostContractionActivity3(i,jj)=(MaxP/Baseline2-1)*100;
     
     Bsl=(Baseline2/Baseline-1)*100;
     Sig=(Fluo/Baseline-1)*100;
     line([-5+locs(jj)+ TimeFlashContra,-3+locs(jj)+ TimeFlashContra], [Max+Bsl Max+Bsl],'Color','k','Linewidth',3);
     line([PosMP+locs(jj)+ TimeFlashContra,PosMP+3+locs(jj)+ TimeFlashContra], [Max+Sig Max+Sig],'Color','r','Linewidth',3);
 end
 
 Max=Max+25;
 
 title (strcat('Cell N°',num2str(i),' / Average Fluo:',num2str(mean(PostContractionActivity(i,:)))))
 
end
        
        title('DF/F versus time for all cells (ROI)')
        xlabel('Time (Frame)')
        ylabel('DF/F ') 
        set(gca,'FontSize',20,'FontWeight','bold')
        set(gcf,'Units','Normalized','Position', [0.05, 0.05, 0.9, 0.85]) 
        saveas(gcf,strcat(folder,image(1:end-4),'_AllSignals.png'));
  
        %Save 
        output.PostContractionActivity=PostContractionActivity;
        output.PostContractionActivity2=PostContractionActivity2;
        output.PostContractionActivity3=PostContractionActivity3;
        output.MeanDFperROI=mean(PostContractionActivity2,2);
        output.SEMDFperROI=std(PostContractionActivity2,0,2)./sqrt(size(PostContractionActivity2,2));
        output.raw_G=raw_G;
        output.dff_G=dff_G;
        output.BaselineROI=BaselineROI;
        output.ref=ref;
        output.TransXY=TransXY;
        output.MotionArtifact=PeaksPos;
        
        output.Gr=Gr;
        output.Go=Go;
        output.Gr=Gr;
        output.rois=rois;
        output.rectindX=rectindX;
        output.rectindY=rectindY;
        output.template=template;
        
        save(strcat(folder,'Analysis_',movieG(1:end-4),'.mat'),'output');
        
end



  