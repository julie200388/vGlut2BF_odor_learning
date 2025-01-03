%% Subtract MO responsive cells
submouseztracesmeanpercell=mouseztracesmeanpercell;submouseztracesmeanpercellaf=mouseztracesmeanpercellaf;
for i=1:12
    for j=1:size(RCall{1},1)
        if RCall{1}(j,1)~=0
           submouseztracesmeanpercell{i}(:,RCall{10}(j,1))=mouseztracesmeanpercell{i}(:,RCall{10}(j,1))-mouseztracesmeanpercell{10}(:,RCall{10}(j,1));
      
        else
        submouseztracesmeanpercell{i}(:,RCall{10}(j,1))=mouseztracesmeanpercell{i}(:,RCall{10}(j,1));
        end
    end
   %  for j=1:size(RCallaf{1},1)
   %      if RCallaf{1}(j,1)~=0
   %         submouseztracesmeanpercellaf{i}(:,RCallaf{10}(j,1))=mouseztracesmeanpercellaf{i}(:,RCallaf{1}(j,1))-mouseztracesmeanpercellaf{1}(:,RCallaf{1}(j,1));
   %      else
   %         submouseztracesmeanpercellaf{i}(:,RCallaf{10}(j,1))=mouseztracesmeanpercellaf{i}(:,RCallaf{1}(j,1));
   %      end
   % 
   % end

end
for i=1:12

        
    for j=1:size(ICall{1},1)
        if ICall{1}(j,1)~=0
           submouseztracesmeanpercell{i}(:,ICall{10}(j,1))=mouseztracesmeanpercell{i}(:,ICall{10}(j,1))-mouseztracesmeanpercell{1}(:,ICall{10}(j,1));
        else
            submouseztracesmeanpercell{i}(:,ICall{10}(j,1))=mouseztracesmeanpercell{i}(:,ICall{10}(j,1));
        end
    end
    % for j=1:size(ICallaf{1},1)
    %     if ICallaf{1}(j,1)~=0
    %        submouseztracesmeanpercellaf{i}(:,ICallaf{1}(j,1))=mouseztracesmeanpercellaf{i}(:,ICallaf{1}(j,1))-mouseztracesmeanpercellaf{1}(:,ICallaf{1}(j,1));
    %     else
    %     submouseztracesmeanpercellaf{i}(:,ICallaf{1}(j,1))=mouseztracesmeanpercellaf{i}(:,ICallaf{1}(j,1));
    %     end
    % end

end
%% plot the traces of all the animals
%for i=1:size(submouseztracesmeanpercell,2)
figure(4) 
plot(Time4plot,mean(submouseztracesmeanpercell{4}(:,142:end),2),'LineWidth',2,'color','k');
hold on
SE = std(submouseztracesmeanpercell{4}(:,142:end)')/sqrt(size(submouseztracesmeanpercell{4}(:,142:end),2)); %I am not sure about what the SE will be??(:,RCall{i})
    CIplus = mean(submouseztracesmeanpercell{4}(:,142:end),2)+(1.96*SE');%(:,RCall{i})
    CIminus = mean(submouseztracesmeanpercell{4}(:,142:end),2)-(1.96*SE');%(:,RCall{i})
  shade(Time4plot,CIplus,Time4plot,CIminus,'FillType',[1 2;2 1],'LineStyle','none','FillColor',[0.8 0.8 0.8])%red shadow [0.6,0,0]
    clear CIplus CIminus SE  
xline(0,'LineWidth',1) %on line
%ax.XAxis.FontSize=30
xline(2,'LineWidth',1) % off line % 1 sec odor presentation
ylim([-0.5 1])
%clear mean
%end
%%
%for i=1:size(submouseztracesmeanpercellaf,2)
figure(4)    
   
plot(Time4plot,mean(submouseztracesmeanpercellaf{4}(:,142:end),2),'LineWidth',2,'color','r');
hold on
SE = std(submouseztracesmeanpercellaf{4}(:,142:end)')/sqrt(size(submouseztracesmeanpercellaf{4}(:,142:end),2)); %I am not sure about what the SE will be??(:,RCall{i})
    CIplus = mean(submouseztracesmeanpercellaf{4}(:,142:end),2)+(1.96*SE');%(:,RCall{i})
    CIminus = mean(submouseztracesmeanpercellaf{4}(:,142:end),2)-(1.96*SE');%(:,RCall{i})
  shade(Time4plot,CIplus,Time4plot,CIminus,'FillType',[1 2;2 1],'LineStyle','none','FillColor',[0.6,0,0])%red shadow [0.6,0,0]
    clear CIplus CIminus SE  
xline(0,'LineWidth',1) %on line
xline(2,'LineWidth',1) % off line % 1 sec odor presentation
%ax.XAxis.FontSize=30
ylim([-0.5 1])
%clear mean
%end

%% plot the heatmap
TF=BS*FR+RS*FR+1;%TF=total frames
mouseztracesmeanHF=mean(submouseztracesmeanpercell{2},1);
%before heatmap
for i=1:size(submouseztracesmeanpercell,2)
    figure(i+10)
    subplot(1,2,1)
    ztracesmeannum=cat(1,submouseztracesmeanpercell{i},mouseztracesmeanHF);
    ztracesmeannum2=sortrows(ztracesmeannum',197,'descend');
    ztracesmeannum=ztracesmeannum2(:,1:196);
h=heatmap(ztracesmeannum,'Colormap',parula,'GridVisible','off');
caxis([-1 2]);
%h.ColorScaling = 'scaledrows';
tickdisplay=strings([TF,1]); %15hz-211/20hz-281
        tickdisplay(BS*FR, 1)='on';
        tickdisplay((BS+2)*FR, 1)='off'; 
tickdisplay2=strings([size(ztracesmeannum,1),1]); 
        tickdisplay2(25, 1)='25';
        tickdisplay2(50, 1)='50';
        tickdisplay2(75, 1)='75';
  tickdisplay2(100, 1)='100'; 
         tickdisplay2(125, 1)='125'; 
         tickdisplay2(150, 1)='150'; 
        tickdisplay2(175, 1)='175'; 
        tickdisplay2(200, 1)='200';
        tickdisplay2(225, 1)='225'; 
        tickdisplay2(250, 1)='250'; 

S=struct(h);
ax=S.Axes;

xline(ax, [BS*FR+0.5],'k','LineWidth',2);
xline(ax, [(BS+2)*FR+0.5],'k','LineWidth',2);
h.XDisplayLabels=tickdisplay;
h.YDisplayLabels=tickdisplay2;
clear h caxis ztracesmeannum ztracesmeannum2
end
%%
%after heatmap
for i=1:size(submouseztracesmeanpercellaf,2)
    figure(i+10)
    subplot(1,2,2)
    ztracesmeannumaf=cat(1,submouseztracesmeanpercellaf{i},mouseztracesmeanHF);
    ztracesmeannum2af=sortrows(ztracesmeannumaf',197,'descend');
    ztracesmeannumaf=ztracesmeannum2af(:,1:196);
h=heatmap(ztracesmeannumaf,'Colormap',parula,'GridVisible','off');
caxis([-1 2]);
%h.ColorScaling = 'scaledrows';
tickdisplay=strings([TF,1]); %15hz-211/20hz-281
        tickdisplay(BS*FR, 1)='on';
        tickdisplay((BS+2)*FR, 1)='off'; 
tickdisplay2=strings([size(ztracesmeannumaf,1),1]); 
        tickdisplay2(25, 1)='25';
        tickdisplay2(50, 1)='50';
        tickdisplay2(75, 1)='75';
        tickdisplay2(100, 1)='100'; 
         tickdisplay2(125, 1)='125'; 
         tickdisplay2(150, 1)='150'; 
        tickdisplay2(175, 1)='175'; 
        tickdisplay2(200, 1)='200'; 
        tickdisplay2(225, 1)='225'; 
        tickdisplay2(250, 1)='250'; 
S=struct(h);
ax=S.Axes;
xline(ax, [BS*FR+0.5],'k','LineWidth',2);
xline(ax, [(BS+2)*FR+0.5],'k','LineWidth',2);
%yline(ax,8.5,'k','LineWidth',1);
%yline(ax,24.5,'k','LineWidth',1);
h.XDisplayLabels=tickdisplay;
h.YDisplayLabels=tickdisplay2;
%clear h caxis ztracesmeannumaf ztracesmeannum2af
end
%% Use three std to find the responsive cells
NO=size(mouseztracesmeanpercell,2);
zAVGbaseline{NO}=[];zAVGResponse{NO}=[];zStdbaseline{NO}=[];AUCabs{NO}=[];
for i=1:NO
    
    zAVGbaseline{i}=mean(submouseztracesmeanpercell{i}(1*FR+1:3*FR,:));%Change the number if you want to define the baseline in a different time frame
    zAVGResponse{i}=mean(submouseztracesmeanpercell{i}(BS*FR+1:(BS+2)*FR,:));
   
    AUCabs{i}=sum(abs(submouseztracesmeanpercell{i}(BS*FR+1:(BS+2)*FR,:)));
    zStdbaseline{i}=std(submouseztracesmeanpercell{i}(1*FR+1:3*FR,:));

    %zAVGResponseall=cat(1,zAVGResponseall,zAVGResponse{i});
    clear n s CN4t 
end
%%
RCall{NO}=[];NCall{NO}=[];ICall{NO}=[];
for i=1:NO
     
     RC=[];NC=[];IC=[];
     for j=1:size(submouseztracesmeanpercell{i},2)
         if zAVGResponse{i}(j)>3*zStdbaseline{i}(j)
             RCstd=j;
             RC=cat(1,RC,RCstd);
             clear RCstd
         elseif zAVGResponse{i}(j)<-3*zStdbaseline{i}(j)
              ICstd=j;   
              IC=cat(1,IC,ICstd);
              clear ICstd
         else 
             NCstd=j;
             NC=cat(1,NC,NCstd);
             clear NCstd
         end
   
     end
     if isempty(RC)==1
         RC=0;
     end
     if isempty(NC)==1
         NC=0;
     end
     if isempty(IC)==1
         IC=0;
     end
    
     RCall{i}=RC;
     NCall{i}=NC;
     ICall{i}=IC;
end
%% Use three std to find the responsive cells after condition/stimulation
NO=size(mouseztracesmeanpercellaf,2);
zAVGbaselineaf{NO}=[];zAVGResponseaf{NO}=[];zStdbaselineaf{NO}=[];zAVGResponseallaf=[];AUCaf{NO}=[];
for i=1:NO
    
    zAVGbaselineaf{i}=mean(submouseztracesmeanpercellaf{i}(1*FR+1:3*FR,:));%Change the number if you want to define the baseline in a different time frame
    zAVGResponseaf{i}=mean(submouseztracesmeanpercellaf{i}(BS*FR+1:(BS+2)*FR,:));

    AUCaf{i}=sum(submouseztracesmeanpercellaf{i}(BS*FR+1:(BS+2)*FR,:));
    zStdbaselineaf{i}=std(submouseztracesmeanpercellaf{i}(1*FR+1:3*FR,:));
    %zAVGResponseallaf=cat(1,zAVGResponseallaf,zAVGResponseaf{i});
    clear n s CN4t 
end
%%
RCallaf{NO}=[];NCallaf{NO}=[];ICallaf{NO}=[];
for i=1:NO
     s=size(submouseztracesmeanpercellaf{i},2);
     RC=[];NC=[];IC=[];
     for j=1:s
         if zAVGResponseaf{i}(j)>3*zStdbaselineaf{i}(j)
             RCstd=j;
             RC=cat(1,RC,RCstd);
             clear RCstd
         elseif zAVGResponseaf{i}(j)<-3*zStdbaselineaf{i}(j)
              ICstd=j;   
              IC=cat(1,IC,ICstd);
              clear ICstd
         else 
             NCstd=j;
             NC=cat(1,NC,NCstd);
             clear NCstd
         end
   
     end
     if isempty(RC)==1
         RC=0;
     end
     if isempty(NC)==1
         NC=0;
     end
     if isempty(IC)==1
         IC=0;
     end
    
     RCallaf{i}=RC;
     NCallaf{i}=NC;
     ICallaf{i}=IC;
end
%% Find RC IC NC in Pentanol and Hexanol
NACmatrixall=[];

    NACmatrix=zeros(NAC,2);
    for i= 2:3
        NACmatrix([RCallaf{i}],i-1)=1;
        NACmatrix([ICallaf{i}],i-1)=-1;
        
    end
    NACmatrix(:,3)=NACmatrix(:,2)-NACmatrix(:,1);
    NACmatrixall=cat(1,NACmatrixall,NACmatrix);
    %%
% In RC to Pentanol

down4FS=size(find(NACmatrixall(:,3)==-2&NACmatrixall(:,1)==1),1);
up4FS=size(find(NACmatrixall(:,3)==0&NACmatrixall(:,1)==1),1);
None4FS=size(find(NACmatrixall(:,3)==-1&NACmatrixall(:,1)==1),1);
RC2Pent=cat(1,down4FS,up4FS,None4FS);
clear down4FS up4FS None4FS

% In IC to Pentanol

down4FS_2=size(find(NACmatrixall(:,3)==0&NACmatrixall(:,1)==-1),1);
up4FS_2=size(find(NACmatrixall(:,3)==2&NACmatrixall(:,1)==-1),1);
None4FS_2=size(find(NACmatrixall(:,3)==1&NACmatrixall(:,1)==-1),1);
IC2Pent=cat(1,down4FS_2,up4FS_2,None4FS_2);
clear down4FS_2 up4FS_2 None4FS_2
% In None to Pentanol

down4FS_3=size(find(NACmatrixall(:,3)==-1&NACmatrixall(:,1)==0),1);
up4FS_3=size(find(NACmatrixall(:,3)==1&NACmatrixall(:,1)==0),1);
None4FS_3=size(find(NACmatrixall(:,3)==0&NACmatrixall(:,1)==0),1);
None2Pent=cat(1,down4FS_3,up4FS_3,None4FS_3);
clear down4FS_3 up4FS_3 None4FS_3
Response2Pent=cat(2,RC2Pent,IC2Pent,None2Pent);
%% Find RC IC NC bf af condition in all animals
NACmatrixallbfaf=[];
for i=2:3
   
    
    NACmatrixallbf=zeros(NAC,1);
    NACmatrixallbf([RCall{i}],1)=1;
    NACmatrixallbf([ICall{i}],1)=-1;

    NACmatrixallaf=zeros(NAC,1);
    NACmatrixallaf([RCallaf{i}],1)=1;
    NACmatrixallaf([ICallaf{i}],1)=-1;
    
    NACmatrixall=cat(2,NACmatrixallbf,NACmatrixallaf);
    NACmatrixallbfaf=cat(2,NACmatrixallbfaf,NACmatrixall);

end

NACmatrixallbfaf(:,5)=NACmatrixallbfaf(:,2)-NACmatrixallbfaf(:,1);
% NACmatrixall([18 22 25 31 45 50 51 58 59 60 65 70 71 76 77 78 79 80 85 88 89 107 119],5)=-3;
% NACmatrixall([28 52 53 97 98 99 100 134 135],5)=-4;
NACmatrixallbfaf(:,6)=NACmatrixallbfaf(:,4)-NACmatrixallbfaf(:,3);
% NACmatrixall([18 22 25 31 45 50 51 58 59 60 65 70 71 76 77 78 79 80 85 88 89 107 119],6)=-3;
% NACmatrixall([28 52 53 97 98 99 100 134 135],6)=-4;
%% 
all{2}=[];
for i=0:1
% In Excited cells before condition

down4FS=size(find(NACmatrixallbfaf(:,5+i)==-2&NACmatrixallbfaf(:,1+2*i)==1),1);
up4FS=size(find(NACmatrixallbfaf(:,5+i)==0&NACmatrixallbfaf(:,1+2*i)==1),1);
None4FS=size(find(NACmatrixallbfaf(:,5+i)==-1&NACmatrixallbfaf(:,1+2*i)==1),1);
%dropout=size(find(NACmatrixall(:,5+i)==-3&NACmatrixall(:,1+2*i)==1),1);

% In Inhibited cells before condition

down4FS_2=size(find(NACmatrixallbfaf(:,5+i)==0&NACmatrixallbfaf(:,1+2*i)==-1),1);
up4FS_2=size(find(NACmatrixallbfaf(:,5+i)==2&NACmatrixallbfaf(:,1+2*i)==-1),1);
None4FS_2=size(find(NACmatrixallbfaf(:,5+i)==1&NACmatrixallbfaf(:,1+2*i)==-1),1);
%dropout_2=size(find(NACmatrixall(:,5+i)==-3&NACmatrixall(:,1+2*i)==-1),1);

% In No Change before condition

down4FS_3=size(find(NACmatrixallbfaf(:,5+i)==-1&NACmatrixallbfaf(:,1+2*i)==0),1);
up4FS_3=size(find(NACmatrixallbfaf(:,5+i)==1&NACmatrixallbfaf(:,1+2*i)==0),1);
None4FS_3=size(find(NACmatrixallbfaf(:,5+i)==0&NACmatrixallbfaf(:,1+2*i)==0),1);
% dropout_3=size(find(NACmatrixall(:,5+i)==-3&NACmatrixall(:,1+2*i)==0),1);
% % Show up
% down_4=size(find(NACmatrixall(:,5+i)==-4&NACmatrixall(:,2+2*i)==-1),1);
% up_4=size(find(NACmatrixall(:,5+i)==-4&NACmatrixall(:,2+2*i)==1),1);
% None_4=size(find(NACmatrixall(:,5+i)==-4&NACmatrixall(:,2+2*i)==0),1);

all{i+1}=cat(1,down4FS,up4FS,None4FS,down4FS_2,up4FS_2,None4FS_2,down4FS_3,up4FS_3,None4FS_3);
end

%% Make the traces ready for PCA
PCAmousetracesmeanpercellbf=[];PCAmousetracesmeanpercellaf=[];
for i=1:5
    PCAmousetracesmeanpercellbf=cat(1,PCAmousetracesmeanpercellbf,submouseztracesmeanpercell{i}(:,104:end));
    PCAmousetracesmeanpercellaf=cat(1,PCAmousetracesmeanpercellaf,submouseztracesmeanpercellaf{i}(:,104:end));

end

%% Before and after condition (Pentanol vs Octanol)
PCAbfaf_pent=cat(1,PCAmousetracesmeanpercellbf(197:392,:),PCAmousetracesmeanpercellbf(785:980,:),PCAmousetracesmeanpercellaf(197:392,:),PCAmousetracesmeanpercellaf(785:980,:));
%Hexbf Octbf Hexaf Octaf
%% Before and after condition (Hexanol vs Octanol)
PCAbfaf_Hex=cat(1,PCAmousetracesmeanpercellbf(393:588,:),PCAmousetracesmeanpercellbf(785:980,:),PCAmousetracesmeanpercellaf(393:588,:),PCAmousetracesmeanpercellaf(785:980,:));
%%
PCAbfaf_pent_Hex_Oct=cat(1,PCAmousetracesmeanpercellbf(1:196,:),PCAmousetracesmeanpercellbf(197:392,:),PCAmousetracesmeanpercellbf(589:784,:),PCAmousetracesmeanpercellaf(1:196,:),PCAmousetracesmeanpercellaf(197:392,:),PCAmousetracesmeanpercellaf(589:784,:));
%% Before and after condition (Pentanol vs Octanol)
CV1 = cov(PCAbfaf_pent);
% compute the eigenvalues and eigenvectors of the covariance matrix
[ev, el] = eig(CV1);
figure;plot(diag(el)./sum(diag(el)),'.')

var = diag(el)./sum(diag(el));
for i = 1:length(var)
    cumvar(i) = sum(var(end-(i-1):end));
end
figure; plot(cumvar)
tmp = find(cumvar > .9);
dimensionality = min(tmp)

% View the eigenvectors across glomeruli (only top 3)
figure;imagesc(ev(:,end-3:end))


% PCA dimensionality reduction
% Project the data to top 3 eigenvectors 
% corresponding to the largest 3 eigenvalues
%data_2d = dataodorsandoffset*ev(:,size(ev)-1:size(ev));
data_3d = PCAbfaf_pent*ev(:,size(ev)-2:size(ev));
data_10d = PCAbfaf_pent*ev(:,size(ev)-9:size(ev));
% Visualize the dimensionality reduced dataset

clr = [0,1,0;0.06,1,1;0,0.3,0.2;0,0,1];%light B light G dark B dark G
% check color scheme: c = uisetcolor
%clr = jet(5);
% Plot the data along the first 3 PCs
% Color the data points based on odor labels
ind=size(PCAbfaf_pent,1)/4;
figure(i+5);hold on
for ii = 1:4
    %ind = find(clab == ii);
%     subplot(3,5,ii);
    plot3(data_3d(1+ind*(ii-1):ind*ii,3),data_3d(1+ind*(ii-1):ind*ii,2),data_3d(1+ind*(ii-1):ind*ii,1),'-','color',clr(ii,:),'linewidth',2,'MarkerEdgeColor',clr(ii,:))
    pause(2)
end
grid on 
xlabel('PC1');

ylabel('PC2');
    zlabel('PC3');
    view([210 30]);

%% Before and after condition (Hexanol vs Octanol)
CV1 = cov(PCAbfaf_Hex);
% compute the eigenvalues and eigenvectors of the covariance matrix
[ev, el] = eig(CV1);
figure;plot(diag(el)./sum(diag(el)),'.')

var = diag(el)./sum(diag(el));
for i = 1:length(var)
    cumvar(i) = sum(var(end-(i-1):end));
end
figure; plot(cumvar)
tmp = find(cumvar > .9);
dimensionality = min(tmp)

% View the eigenvectors across glomeruli (only top 3)
figure;imagesc(ev(:,end-3:end))


% PCA dimensionality reduction
% Project the data to top 3 eigenvectors 
% corresponding to the largest 3 eigenvalues
%data_2d = dataodorsandoffset*ev(:,size(ev)-1:size(ev));
data_3d = PCAbfaf_Hex*ev(:,size(ev)-2:size(ev));
data_10d = PCAbfaf_Hex*ev(:,size(ev)-9:size(ev));
% Visualize the dimensionality reduced dataset
clr = [0.92,0.61,0.61;0.06,1,1;0.64,0.08,0.18;0,0,1];
%clr = [0.06,1,1;0,1,0;0,0,1;0,0.3,0.2];
% check color scheme: c = uisetcolor
%clr = jet(5);
% Plot the data along the first 3 PCs
% Color the data points based on odor labels
ind=size(PCAbfaf_Hex,1)/4;
figure;hold on
for ii = 1:4
    %ind = find(clab == ii);
%     subplot(3,5,ii);
    plot3(data_3d(1+ind*(ii-1):ind*ii,3),data_3d(1+ind*(ii-1):ind*ii,2),data_3d(1+ind*(ii-1):ind*ii,1),'-','color',clr(ii,:),'linewidth',2,'MarkerEdgeColor',clr(ii,:))
    pause(2)
end
grid on 
xlabel('PC1');

ylabel('PC2');
    zlabel('PC3');
    view([210 30]);
    %% Before and after condition (Pentanol vs Octanol)
CV1 = cov(PCAbfaf_pent_Hex_Oct);
% compute the eigenvalues and eigenvectors of the covariance matrix
[ev, el] = eig(CV1);
figure;plot(diag(el)./sum(diag(el)),'.')

var = diag(el)./sum(diag(el));
for i = 1:length(var)
    cumvar(i) = sum(var(end-(i-1):end));
end
figure; plot(cumvar)
tmp = find(cumvar > .9);
dimensionality = min(tmp)

% View the eigenvectors across glomeruli (only top 3)
figure;imagesc(ev(:,end-3:end))


% PCA dimensionality reduction
% Project the data to top 3 eigenvectors 
% corresponding to the largest 3 eigenvalues
%data_2d = dataodorsandoffset*ev(:,size(ev)-1:size(ev));
data_3d = PCAbfaf_pent_Hex_Oct*ev(:,size(ev)-2:size(ev));
data_10d = PCAbfaf_pent_Hex_Oct*ev(:,size(ev)-9:size(ev));
% Visualize the dimensionality reduced dataset

clr = [0,1,0;0.92,0.61,0.61;0.06,1,1;0,0.3,0.2;0.64,0.08,0.18;0,0,1];%light B light G dark B dark G
% check color scheme: c = uisetcolor
%clr = jet(5);
% Plot the data along the first 3 PCs
% Color the data points based on odor labels
ind=size(PCAbfaf_pent_Hex_Oct,1)/6;
figure();hold on
for ii = 1:6
    % if ii==2
    %     continue
    % end
    % if ii==5 
    %     continue
    % end
   % if ii==4 
   %      continue
   %  end
    %ind = find(clab == ii);
%     subplot(3,5,ii);
    plot3(data_3d(1+ind*(ii-1):ind*ii,3),data_3d(1+ind*(ii-1):ind*ii,2),data_3d(1+ind*(ii-1):ind*ii,1),'-','color',clr(ii,:),'linewidth',2,'MarkerEdgeColor',clr(ii,:))
    pause(2)
end
grid on 
xlabel('PC1');
xlim([-2 14]);
ylabel('PC2');
ylim([-4 8]);
    zlabel('PC3');
   zlim([-6 6]); 
    view([210 30]);

    %%
    PCAmousetracesmeanpercellbf=PCAmousetracesmeanpercellbf(197:end,:);
    PCAmousetracesmeanpercellaf=PCAmousetracesmeanpercellaf(197:end,:);
    %% Before and after condition (Pentanol vs Octanol)
CV1 = cov(PCAmousetracesmeanpercellbf);
% compute the eigenvalues and eigenvectors of the covariance matrix
[ev, el] = eig(CV1);
figure(1);plot(diag(el)./sum(diag(el)),'.')

var = diag(el)./sum(diag(el));
for i = 1:length(var)
    cumvar(i) = sum(var(end-(i-1):end));
end
figure(2); plot(cumvar)
tmp = find(cumvar > .9);
dimensionality = min(tmp)

% View the eigenvectors across glomeruli (only top 3)
figure(3);imagesc(ev(:,end-3:end))


% PCA dimensionality reduction
% Project the data to top 3 eigenvectors 
% corresponding to the largest 3 eigenvalues
%data_2d = dataodorsandoffset*ev(:,size(ev)-1:size(ev));
data_3d = PCAmousetracesmeanpercellbf*ev(:,size(ev)-2:size(ev));
data_10d = PCAmousetracesmeanpercellbf*ev(:,size(ev)-9:size(ev));
% Visualize the dimensionality reduced dataset
clr = [0, 1, 0;   % Green
       1, 0, 0;   % Red
       1, 1, 0;   % Yellow
       0, 0, 1];  % Blue
%clr = [0,1,0;0.06,1,1;0,0.3,0.2;0,0,1];%light B light G dark B dark G
% check color scheme: c = uisetcolor
%clr = jet(5);
% Plot the data along the first 3 PCs
% Color the data points based on odor labels
ind=size(PCAmousetracesmeanpercellbf,1)/4;
figure(4);hold on
for ii = 1:4
    %ind = find(clab == ii);
%     subplot(3,5,ii);
    plot3(data_3d(1+ind*(ii-1):ind*ii,3),data_3d(1+ind*(ii-1):ind*ii,2),data_3d(1+ind*(ii-1):ind*ii,1),'-','color',clr(ii,:),'linewidth',2,'MarkerEdgeColor',clr(ii,:))
    pause(2)
end
grid on 
xlabel('PC1');
xlim([-14 2]); % Corrected to ascending order
ylabel('PC2');
ylim([-4 8]); % Corrected to ascending order
zlabel('PC3');
zlim([-4 8]);
view([210 30]);
   %% Before and after condition (Pentanol vs Octanol)
CV1 = cov(PCAmousetracesmeanpercellaf);
% compute the eigenvalues and eigenvectors of the covariance matrix
[ev, el] = eig(CV1);
figure(5);plot(diag(el)./sum(diag(el)),'.')

var = diag(el)./sum(diag(el));
for i = 1:length(var)
    cumvar(i) = sum(var(end-(i-1):end));
end
figure(6); plot(cumvar)
tmp = find(cumvar > .9);
dimensionality = min(tmp)

% View the eigenvectors across glomeruli (only top 3)
figure(7);imagesc(ev(:,end-3:end))


% PCA dimensionality reduction
% Project the data to top 3 eigenvectors 
% corresponding to the largest 3 eigenvalues
%data_2d = dataodorsandoffset*ev(:,size(ev)-1:size(ev));
data_3d = PCAmousetracesmeanpercellaf*ev(:,size(ev)-2:size(ev));
data_10d = PCAmousetracesmeanpercellaf*ev(:,size(ev)-9:size(ev));
% Visualize the dimensionality reduced dataset

%clr = [0,1,0;0.06,1,1;0,0.3,0.2;0,0,1];%light B light G dark B dark G
% check color scheme: c = uisetcolor
clr = [0, 1, 0;   % Green
       1, 0, 0;   % Red
       1, 1, 0;   % Yellow
       0, 0, 1];  % Blue
% Plot the data along the first 3 PCs
% Color the data points based on odor labels
ind=size(PCAmousetracesmeanpercellaf,1)/4;
figure(8);hold on
for ii = 1:4
    %ind = find(clab == ii);
%     subplot(3,5,ii);
    plot3(data_3d(1+ind*(ii-1):ind*ii,3),data_3d(1+ind*(ii-1):ind*ii,2),data_3d(1+ind*(ii-1):ind*ii,1),'-','color',clr(ii,:),'linewidth',2,'MarkerEdgeColor',clr(ii,:))
    pause(2)
end
grid on 
xlabel('PC1');
ylabel('PC2');
    zlabel('PC3');
    view([210 30]);



 %% Calculate the angle between before and after condition
mouseztraces{5}=[];mouseztracesaf{5}=[];

for i=1:5
    for mouse=1:8
        mouseztraces{i}=cat(2,mouseztraces{i},mousenum{1,mouse}.ztemptraceall{i});
        mouseztracesaf{i}=cat(2,mouseztracesaf{i},mousenum{1,mouse}.ztemptraceallaf{i});
    end
end
%%
samplingRate = 15;  % Frames per second, adjust based on your data
duration = 2;  % Duration in seconds for the odor presentation window
timeWindow = 76:90;  % Time points for the first 2 seconds
% Initialize array to store angles for each trial
anglestoControl = zeros(154, 4);
distancetoControl=zeros(154, 4);
% Loop through each trial to calculate the angle
for i =1:4
for neuron = 1:154
    % Extract calcium activity for the first 2 seconds (time window)
    activityWindow = mouseztraces{i}(timeWindow, neuron, :);
    activityWindowaf=mouseztraces{5}(timeWindow, neuron, :);
    % Calculate the mean activity vector over the time window
    meanActivity = mean(activityWindow, 3);  % Mean across time points
    meanActivityaf = mean(activityWindowaf, 3); 
 
    % Calculate the dot product between the two mean activity vectors
    dotProduct = dot(meanActivity, meanActivityaf);
    % ------ Euclidean Distance ------
        euclideanDistance = norm(meanActivityaf - meanActivity);  % Euclidean distance
    % Calculate the norms (magnitudes) of both activity vectors
    normActivity1 = norm(meanActivity);  % Magnitude of the first activity vector
    normActivity2 = norm(meanActivityaf);  % Magnitude of the second activity vector
    
    % Calculate the angle in radians between the two vectors using acos
    angle = acos(dotProduct / (normActivity1 * normActivity2));  % Angle in radians
    
    % Store the angle for this trial
    anglestoControl (neuron,i) = angle*180/pi;
    distancetoControl(neuron,i)=euclideanDistance;
end
end
%%
samplingRate = 1;  % Frames per second, adjust based on your data
windowDuration = 1;  % Duration in seconds for the odor presentation window
windowSize = windowDuration * samplingRate;  % Calculate the number of frames in the time window
numTimePoints = 196;  % Total number of time points in the data (adjust based on your data)
numComparisons = numTimePoints - windowSize;  % Number of consecutive time window comparisons

% Initialize arrays to store angles and distances for each neuron, trial, and comparison
anglesBetweenSets = zeros(154, 900, 5);
distanceBetweenSets = zeros(154, 900, 5);

% Loop through each trial to calculate the angle and distance
for i = 1:5
    for j =1:5
    for neuron = 1:154
        % Loop through consecutive time windows for each trial
        for t = 1:30  % Loop from the first to the penultimate window
            % Extract calcium activity for the current and previous time windows
            currentWindow = mouseztraces{i}(t+44:(t+44+windowSize), neuron, :);
            previousWindow = mouseztraces{j}(t+44+30:(t+44+30+windowSize), neuron, :);
            
            % Calculate the mean activity vector over the time window
            meanCurrentActivity = mean(currentWindow, 3);  % Mean across time points
            meanPreviousActivity = mean(previousWindow, 3);  
            
            % Calculate the dot product between the two mean activity vectors
            dotProduct = dot(meanCurrentActivity, meanPreviousActivity);
            
            % ------ Euclidean Distance ------
            euclideanDistance = norm(meanPreviousActivity - meanCurrentActivity);  % Euclidean distance
            
            % Calculate the norms (magnitudes) of both activity vectors
            normCurrentActivity = norm(meanCurrentActivity);  % Magnitude of the current activity vector
            normPreviousActivity = norm(meanPreviousActivity);  % Magnitude of the previous activity vector
            
            % Calculate the angle in radians between the two vectors using acos
            angle = acos(dotProduct / (normCurrentActivity * normPreviousActivity));  % Angle in radians
            
            % Store the angle and Euclidean distance for this comparison
            anglesBetweenSets(neuron, 30*(t-1)+q, j) = angle * 180 / pi;  % Convert to degrees
            distanceBetweenSets(neuron, 30*(t-1)+q,j) = euclideanDistance;
        end
       
        
    end
end
%%
samplingRate = 1;  % Frames per second, adjust based on your data
windowDuration = 1;  % Duration in seconds for the odor presentation window
windowSize = windowDuration * samplingRate;  % Calculate the number of frames in the time window
numTimePoints = 196;  % Total number of time points in the data (adjust based on your data)
numComparisons = numTimePoints - windowSize;  % Number of consecutive time window comparisons

% Initialize arrays to store angles and distances for each neuron, trial, and comparison
anglesBetweenSets = zeros(154, 900, 5);
distanceBetweenSets = zeros(154, 900, 5);
medall=zeros(5,5);
% Loop through each trial to calculate the angle and distance
for i = 1:5
    for neuron = 1:154
        % Loop through consecutive time windows for each trial
        for t = 1:30  % Loop from the first to the penultimate window
            % Extract calcium activity for the current and previous time windows
            currentWindow = mouseztraces{i}(t+44:(t+44+windowSize), neuron, :);
            for q=1:30
            previousWindow = mouseztraces{i}(q+44+30:(q+44+30+windowSize), neuron, :);
            
            % Calculate the mean activity vector over the time window
            meanCurrentActivity = mean(currentWindow, 3);  % Mean across time points
            meanPreviousActivity = mean(previousWindow, 3);  
            
            % Calculate the dot product between the two mean activity vectors
            dotProduct = dot(meanCurrentActivity, meanPreviousActivity);
            
            % ------ Euclidean Distance ------
            euclideanDistance = norm(meanPreviousActivity - meanCurrentActivity);  % Euclidean distance
            
            % Calculate the norms (magnitudes) of both activity vectors
            normCurrentActivity = norm(meanCurrentActivity);  % Magnitude of the current activity vector
            normPreviousActivity = norm(meanPreviousActivity);  % Magnitude of the previous activity vector
            
            % Calculate the angle in radians between the two vectors using acos
            angle = acos(dotProduct / (normCurrentActivity * normPreviousActivity));  % Angle in radians
            
            % Store the angle and Euclidean distance for this comparison
            anglesBetweenSets(neuron, 30*(t-1)+q, i) = angle * 180 / pi;  % Convert to degrees
            distanceBetweenSets(neuron, 30*(t-1)+q,i) = euclideanDistance;
            end
        
        end
    end
    medall(i,i)=median(median(distanceBetweenSets(:,:,i)));
    
end
%%
samplingRate = 1;  % Frames per second, adjust based on your data
windowDuration = 1;  % Duration in seconds for the odor presentation window
windowSize = windowDuration * samplingRate;  % Calculate the number of frames in the time window
numTimePoints = 196;  % Total number of time points in the data (adjust based on your data)
numComparisons = numTimePoints - windowSize;  % Number of consecutive time window comparisons

% Initialize arrays to store angles and distances for each neuron, trial, and comparison
anglesBetweenSetsaf = zeros(154, 900, 5);
distanceBetweenSetsaf = zeros(154, 900, 5);
medallaf=zeros(5,900);
% Loop through each trial to calculate the angle and distance
for i = 1:5
    for neuron = 1:154
        % Loop through consecutive time windows for each trial
        for t = 1:30  % Loop from the first to the penultimate window
            % Extract calcium activity for the current and previous time windows
            currentWindow = mouseztracesaf{i}(t+44:(t+44+windowSize), neuron, :);
            for q=1:30
            previousWindow = mouseztracesaf{i}(q+44+30:(q+44+30+windowSize), neuron, :);
            
            % Calculate the mean activity vector over the time window
            meanCurrentActivity = mean(currentWindow, 3);  % Mean across time points
            meanPreviousActivity = mean(previousWindow, 3);  
            
            % Calculate the dot product between the two mean activity vectors
            dotProduct = dot(meanCurrentActivity, meanPreviousActivity);
            
            % ------ Euclidean Distance ------
            euclideanDistance = norm(meanPreviousActivity - meanCurrentActivity);  % Euclidean distance
            
            % Calculate the norms (magnitudes) of both activity vectors
            normCurrentActivity = norm(meanCurrentActivity);  % Magnitude of the current activity vector
            normPreviousActivity = norm(meanPreviousActivity);  % Magnitude of the previous activity vector
            
            % Calculate the angle in radians between the two vectors using acos
            angle = acos(dotProduct / (normCurrentActivity * normPreviousActivity));  % Angle in radians
            
            % Store the angle and Euclidean distance for this comparison
            anglesBetweenSetsaf(neuron, 30*(t-1)+q, i) = angle * 180 / pi;  % Convert to degrees
            distanceBetweenSetsaf(neuron, 30*(t-1)+q,i) = euclideanDistance;
            end
        
        end
    end
    medallaf(i,:)=median(distanceBetweenSetsaf(:,:,i));

end
%%  Coefficient by timepoints_after conditioning/stimulation
Coallmice=[];
for mouse=1:8
Coall=[];cellmean{5}=[];
for i= 1:5
    Cocells=[];
%     for j=1:10
%         if j==1
%         alltrials=idzodortraces{i}(:,1:NAC);
%         else
%         alltrials(:,:,j)=idzodortraces{i}(:,NAC*(j-1)+1:NAC*j);
%         end
%     end
    for j=1:10
        cellmean{i}=mean(mousenum{1,mouse}.ztemptraceall{i}(46:75,:,:),2);%(46:90,:,:)
        cellmean{i}(:,:,11)=cellmean{i}(:,:,j);
        cellmean{i}(:,:,j)=[];
        Co=corrcoef(cellmean{i}(:,:,10),mean(cellmean{i}(:,:,1:9),3));
        Cocells=cat(2,Cocells,Co(2));
        clear Co
    end
 Coall=cat(1,Coall,Cocells);

end
Comean=mean(Coall,2);
Coall=cat(2,Coall,Comean);
Coallmice=cat(1,Coallmice,Coall);
end
%%  Coefficient by timepoints_after conditioning/stimulation
Coallmiceaf=[];
for mouse=1:8
Coallaf=[];cellmeanaf{5}=[];
for i= 1:k
    Cocellsaf=[];
%     for j=1:10
%         if j==1
%         alltrials=idzodortraces{i}(:,1:NAC);
%         else
%         alltrials(:,:,j)=idzodortraces{i}(:,NAC*(j-1)+1:NAC*j);
%         end
%     end
    for j=1:10
        cellmeanaf{i}=mean(mousenum{1,mouse}.ztemptraceallaf{i}(46:75,:,:),2);%(46:90,:,:)
        cellmeanaf{i}(:,:,11)=cellmeanaf{i}(:,:,j);
        cellmeanaf{i}(:,:,j)=[];
        Co=corrcoef(cellmeanaf{i}(:,:,10),mean(cellmeanaf{i}(:,:,1:9),3));
        Cocellsaf=cat(2,Cocellsaf,Co(2));
        clear Co
    end
 Coallaf=cat(1,Coallaf,Cocellsaf);

end
Comeanaf=mean(Coallaf,2);
Coallaf=cat(2,Coallaf,Comeanaf);
Coallmiceaf=cat(1,Coallmiceaf,Coallaf);
end
%% Save figures
%% Save the files for four odors
saveas(figure(1),'12animals_Mineral oil_bafcondition_MOsubtracted.png')
saveas(figure(2),'12animals_Pentanol_bafcondition_MOsubtracte.png')
saveas(figure(3),'12animals_Hexanol_bafcondition_MOsubtracte.png')
saveas(figure(4),'12animals_Heptanol_bafcondition_MOsubtracte.png')
saveas(figure(5),'12animals_Octanol_bafcondition_MOsubtracte.png')
%%
saveas(figure(11),'12animals_Mineral oil_bafcondition_MOsubtracted_heatmap.png')
saveas(figure(12),'12animals_Pentanol_bafcondition_MOsubtracted_heatmap.png')
saveas(figure(13),'12animals_Hexanol_bafcondition_MOsubtracted_heatmap.png')
saveas(figure(14),'12animals_Heptanol_bafcondition_MOsubtracted_heatmap.png')
saveas(figure(15),'12animals_Octanol_bafcondition_MOsubtracted_heatmap.png')
%% Average response of RC, IC, NC
zAVGResponseRC=zeros(2,257,5);
zAVGResponseIC=zeros(2,257,5);
zAVGResponseNC=zeros(2,257,5);
for i=1:5
zAVGResponseRC(1,1:size(RCall{i},1),i)=mean(submouseztracesmeanpercell{i}(BS*FR+1:(BS+2)*FR,RCall{i}));
zAVGResponseIC(1,1:size(ICall{i},1),i)=mean(submouseztracesmeanpercell{i}(BS*FR+1:(BS+2)*FR,ICall{i}));
zAVGResponseNC(1,1:size(NCall{i},1),i)=mean(submouseztracesmeanpercell{i}(BS*FR+1:(BS+2)*FR,NCall{i}));
zAVGResponseRC(2,1:size(RCallaf{i},1),i)=mean(submouseztracesmeanpercellaf{i}(BS*FR+1:(BS+2)*FR,RCallaf{i}));
zAVGResponseIC(2,1:size(ICallaf{i},1),i)=mean(submouseztracesmeanpercellaf{i}(BS*FR+1:(BS+2)*FR,ICallaf{i}));
zAVGResponseNC(2,1:size(NCallaf{i},1),i)=mean(submouseztracesmeanpercellaf{i}(BS*FR+1:(BS+2)*FR,NCallaf{i}));
end
%%
Common_RC_af=intersect(RCallaf{2},RCallaf{3});
Common_IC_af=intersect(ICallaf{2},ICallaf{3});
zAVG_CommonRC=cat(1,zAVGResponseaf{2}(:,Common_RC_af), zAVGResponseaf{3}(:,Common_RC_af));
zAVG_CommonIC=cat(1,zAVGResponseaf{2}(:,Common_IC_af), zAVGResponseaf{3}(:,Common_IC_af));
%%

numSplits = 10;  % Define the number of splits you want
corrMatrix = zeros(4); 
corrMatrixmice=zeros(4,4,8);
for mouse=1:8
for i=2:5
     trialsA=mousenum2{1,mouse}.ztemptraceall{i}(46:end,:,:);
    for t=2:5 
        trialsB=mousenum2{1,mouse}.ztemptraceall{t}(46:end,:,:);
        
        correlations = zeros(1, numSplits);
        for split = 1:numSplits
            % Perform a random split for each iteration
           
            
            numTrialsA = size(trialsA, 3);
            numTrialsB = size(trialsB, 3);
            idxA = randperm(numTrialsA);
            idxB = randperm(numTrialsB);
            splitA1 = mean(trialsA(:, :, idxA(1:floor(numTrialsA/2))), 3); 
            splitB1 = mean(trialsB(:, :, idxB(1:floor(numTrialsB/2))), 3);
            
            vectorA1 = splitA1(:);
            vectorB1 = splitB1(:);
            
            % Calculate correlations for this split
            correlations(split) = corr(vectorA1, vectorB1);
        end
         corrMatrix(i-1,t-1) = mean(correlations);
    end

end
corrMatrixmice(:,:,mouse)=corrMatrix;

end
corrMatrixmicemean=mean(corrMatrixmice,3);
clear corrMatrix correlations
%%
% Plot heat map
figure(1)
imagesc(corrMatrixmicemean);  % Display the correlation matrix as a heat map
colorbar;              % Add a color bar to indicate the correlation scale
title('Odor Correlation Matrix');
xlabel('Odor');
ylabel('Odor');
caxis([0.3 0.9]);% Adjust axis ticks and labels if you have specific odor labels
numOdors = size(corrMatrix, 1);
xticks(1:numOdors);
yticks(1:numOdors);
% Optionally, label the ticks if odors have specific names
odorLabels = {'Pentanol(HF)', 'Hexanol(HF)', 'Heptanol(Ctrl)', 'Octanol(Ctrl)'};
xticklabels(odorLabels);
yticklabels(odorLabels);
for i = 1:4
    for j = 1:4
        % Display the correlation value in each cell
        text(j, i, num2str(corrMatrixmicemean(i, j), '%.2f'), ... % Format to 2 decimal places
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', ...
             'Color', 'white'); % Change color if needed for visibility
    end
end
%%
numSplits = 10;  % Define the number of splits you want
corrMatrix = zeros(4); 
corrMatrixmiceaf=zeros(4,4,8);
for mouse=1:8
for i=2:5
     trialsA=mousenum2{1,mouse}.ztemptraceallaf{i}(46:end,:,:);
    for t=2:5 
        trialsB=mousenum2{1,mouse}.ztemptraceallaf{t}(46:end,:,:);
        
        correlations = zeros(1, numSplits);
        for split = 1:numSplits
            % Perform a random split for each iteration
           
            
            numTrialsA = size(trialsA, 3);
            numTrialsB = size(trialsB, 3);
            idxA = randperm(numTrialsA);
            idxB = randperm(numTrialsB);
            splitA1 = mean(trialsA(:, :, idxA(1:floor(numTrialsA/2))), 3); 
            splitB1 = mean(trialsB(:, :, idxB(1:floor(numTrialsB/2))), 3);
            
            vectorA1 = splitA1(:);
            vectorB1 = splitB1(:);
            
            % Calculate correlations for this split
            correlations(split) = corr(vectorA1, vectorB1);
        end
         corrMatrix(i-1,t-1) = mean(correlations);
    end

end
corrMatrixmiceaf(:,:,mouse)=corrMatrix;

end
corrMatrixmicemeanaf=mean(corrMatrixmiceaf,3);

%%
% Plot heat map
figure(2);
imagesc(corrMatrixmicemeanaf);  % Display the correlation matrix as a heat map
colorbar;              % Add a color bar to indicate the correlation scale
title('Odor Correlation Matrix');
xlabel('Odor');
ylabel('Odor');
caxis([0.3 0.9]);
numOdors = size(corrMatrix, 1);
xticks(1:numOdors);
yticks(1:numOdors);
% Optionally, label the ticks if odors have specific names
odorLabels = {'Pentanol(HF)', 'Hexanol(HF)', 'Heptanol(Ctrl)', 'Octanol(Ctrl)'};
xticklabels(odorLabels);
yticklabels(odorLabels);

for i = 1:4
    for j = 1:4
        % Display the correlation value in each cell
        text(j, i, num2str(corrMatrixmicemeanaf(i, j), '%.2f'), ... % Format to 2 decimal places
             'HorizontalAlignment', 'center', ...
             'VerticalAlignment', 'middle', ...
             'Color', 'white'); % Change color if needed for visibility
    end
end
%%
saveas(figure(1),'8animals_PCA_4odorsbfcondition_1.png')
saveas(figure(2),'8animals_PCA_4odorsbfcondition_2.png')
saveas(figure(3),'8animals_PCA_4odorsbfcondition_3.png')
saveas(figure(5),'8animals_PCA_4odorsafcondition_1.png')
saveas(figure(6),'8animals_PCA_4odorsafcondition_2.png')
saveas(figure(7),'8animals_PCA_4odorsafconditiont_3.png')
