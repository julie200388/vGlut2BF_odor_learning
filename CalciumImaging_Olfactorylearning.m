%% Calcium imaging data processing for individual mouse from the associative learning experiments
%Open calcium imaging trace files
clear all
close all
NAC = input('number of accepted cells: '); %Number of accepted cells
[phf,php] = uigetfile('/Documents/*.csv');
photfn = sprintf('%s/%s',php,phf);
tbl = readtable(photfn,'HeaderLines',0,'TextType','string');
VarNames = tbl.Properties.VariableNames;

allAC=tbl.("accepted");

for i=1:NAC-1 %the number of accepted cells

    AC=tbl.("accepted_"+num2str(i));
    allAC=cat(2,allAC,AC);
    clear AC
end
Time=tbl.("Time_s__CellStatus");

%% Seperate allAC into different days
D2start=find(diff(Time)<-1);
%% Analyzing bf(Day1)
allACbf=allAC(1:(D2start),:);
Timebf=Time(1:(D2start),:);
allACbf(isnan(allACbf)==1)=0;

%% Analyzing Af(Day2)
allACaf=allAC(D2start+1:end,:);
Timeaf=Time(D2start+1:end,:);
allACaf(isnan(allACaf)==1)=0;
%% only analyzing one day_RUN THIS AFTER getting allAC
   Timebf=Time;
   allACbf=allAC;
   allACbf(isnan(allACbf))=0;
%% settining up the DIO files_before
DIO2 = readtable('A138Nt_0821bfcondition_GPIO.csv');
%%
DIO2=table2array(DIO2);
a=find(DIO2(:,2)>1000);
b=find(DIO2(:,2)<1000);
DIO2(a(1:end),2)=1;
DIO2(b(1:end),2)=0;
d=diff(DIO2(:,2));
d2=find(d==1)+1;
odorpresentationtime=DIO2(d2,1);
%% setting up the DIO files_after
DIO3= readtable('A138Nt_0829afcondition_GPIO.csv');
%%
DIO3=table2array(DIO3);
a=find(DIO3(:,2)>1000);
b=find(DIO3(:,2)<1000);
DIO3(a(1:end),2)=1;
DIO3(b(1:end),2)=0;
d=diff(DIO3(:,2));
d2=find(d==1)+1;
odorpresentationtime2=DIO3(d2,1);

%% Combine odor presentation and traces % 011422 data don't run this section
[fn3,pth3] = uigetfile('*.*');
odorsfn = sprintf('%s%s',pth3,fn3);
odors = csvread(odorsfn);
%% %% Combine odor presentation and traces % 011422 data don't run this section
[fn5,pth5] = uigetfile('*.*');
odorsfnaf = sprintf('%s%s',pth5,fn5);
odorsaf = csvread(odorsfnaf);

%% average traces before condition
FR = input('Frame rate: '); %Number of accepted cells
BS=  input('Time before stimulation: ');
RS=  input('Time after stimulation:');
Timebf=Time;allACbf=allAC;
ztemptraceall{max(odors)}=[];
for i = 1:length(odorpresentationtime) %for 011422 file use this:i=1:(length(TP)-1)
    ch = odors(i);
  
    [val,idx] = min(abs(Timebf-odorpresentationtime(i,1)));%This finds the closest actual time value in the trace times that matches the odor time.
    temptrace = allACbf(idx-BS*FR:idx+RS*FR,:);
    activetraces=temptrace-mean(temptrace(FR+1:BS*FR,:),1);
    ztemptrace=activetraces./std(mean(temptrace(FR+1:BS*FR,:),1));%(activetraces-mean(activetraces,1))
    temptracetime = Timebf(idx-BS*FR:idx+RS*FR);
    ztemptraceall{ch}=cat(3,ztemptraceall{ch},ztemptrace);
    clear ch startsec val idx temptrace temptracetime activetraces
end

ztracesmeanpercell{max(odors)}=[];
for i=1:max(odors)
    ztracesmeanpercell{i}=mean(ztemptraceall{i},3);
end
clear ztemptrace


%% average traces after condition
ztemptraceallaf{max(odorsaf)}=[];
for i = 1:length(odorpresentationtime2) %for 011422 file use this:i=1:(length(TP)-1)

    ch = odorsaf(i);
  
    [val,idx] = min(abs(Timeaf-odorpresentationtime2(i,1)));%This finds the closest actual time value in the trace times that matches the odor time.
   
    temptrace = allACaf(idx-BS*FR:idx+RS*FR,:);
    activetraces=temptrace-mean(temptrace(FR+1:BS*FR,:),1);
    ztemptrace=activetraces./std(mean(temptrace(FR+1:BS*FR,:),1));%(activetraces-mean(activetraces,1))
    temptracetimeaf = Timeaf(idx-BS*FR:idx+RS*FR);
    
    ztemptraceallaf{ch}=cat(3,ztemptraceallaf{ch},ztemptrace);

    clear ch startsec val idx temptrace temptracetime activetraces
end

ztracesmeanpercellaf{max(odorsaf)}=[];
for i=1:max(odorsaf)
    ztracesmeanpercellaf{i}=mean(ztemptraceallaf{i},3);
end
%% Output traces 2 seconds after odor onset for the linear model (One file for each mouse)
traces4model2secbf=[];
NO=size(ztemptraceall,2);
for i=1:NO
    featureall=[];featureallaf=[];
for j= 1:10
    feature=reshape(ztemptraceall{i}(46:75,:,j),1,NAC*45);
    featureaf=reshape(ztemptraceallaf{i}(46:75,:,j),1,NAC*45);
    featureall=cat(1,featureall,feature);
    featureallaf=cat(1,featureallaf,featureaf);
end
traces4model2secbf=cat(1,traces4model2secbf,featureall);
traces4model2secaf=cat(1,traces4model2secaf,featureallaf);
odor4feature(10*(i-1)+1:10*i,1)=i;
end

traces4model2secbf=cat(2,traces4model2secbf,odor4feature);
traces4model2secaf=cat(2,traces4model2secaf,odor4feature);
%writematrix(traces4model2secbf,'072624_A134R1_tracesfeature_odorbfcondition_2SECAFODORONSET_v3.csv')
%writematrix(traces4model2secaf,'080324_A134R1_tracesfeature_odorafcondition_2SECAFODORONSET_v3.csv')

% save individual mouse before moving on
%% Open all the files from different animals
% Save the data of individual mouse before this processing
n = input('number of mice: ');
k= input('number of odors: ');
for i = 1:n
    tracevar = uigetfile('*.*');
    mousenum{i} = load(tracevar);
    name1 = strsplit(tracevar, '_');
    name1  = char(name1{1,1});
    mouseID{i} = name1(1:5);
end
%% Combine all the animals_ before condition
mouseztracesmeanpercell{5}=[];
for i=1:5
    for mouse=1:n
        if n==1
       
        
        mouseztracesmeanpercell{i}=mean(mousenum{1,mouse}.ztemptraceall{i},3);
        end
        
        if n>1
      
        mouseztracesmeanpercell{i}=cat(2,mouseztracesmeanpercell{i},mean(mousenum{1,mouse}.ztemptraceall{i},3));  
        end
    end
end
%% Use three std to find the MINERAL OIL responsive cells
NO=size(mouseztracesmeanpercell,2);
zAVGbaseline{NO}=[];zAVGResponse{NO}=[];zStdbaseline{NO}=[];AUCabs{NO}=[];  
    zAVGbaseline{1}=mean(mouseztracesmeanpercell{i}(1*FR+1:3*FR,:));%Change the number if you want to define the baseline in a different time frame
    zAVGResponse{1}=mean(mouseztracesmeanpercell{i}(BS*FR+1:(BS+2)*FR,:));
    zStdbaseline{1}=std(mouseztracesmeanpercell{i}(1*FR+1:3*FR,:));

%%
RCall{NO}=[];NCall{NO}=[];ICall{NO}=[];

     
     RC=[];NC=[];IC=[];
     for j=1:size(mouseztracesmeanpercell{i},2)
         if zAVGResponse{1}(j)>3*zStdbaseline{1}(j)
             RCstd=j;
             RC=cat(1,RC,RCstd);
             clear RCstd
         elseif zAVGResponse{1}(j)<-3*zStdbaseline{1}(j)
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
    
     RCall{1}=RC;
     NCall{1}=NC;
     ICall{1}=IC;
%%
mouseztracesmeanpercellaf{5}=[];
for i=1:5
    for mouse=1:n
        if n==1
        mouseztracesmeanpercellaf{i}=mousenum{1,mouse}.ztracesmeanpercellaf{i};
        end
        
        if n>1
    
        mouseztracesmeanpercellaf{i}=cat(2,mouseztracesmeanpercellaf{i},mousenum{1,mouse}.ztracesmeanpercellaf{i});  
        end
    end
end
%% Use three std to find the MINERAL OIL responsive cells after condition
NO=size(mouseztracesmeanpercellaf,2);
zAVGbaselineaf{NO}=[];zAVGResponseaf{NO}=[];zStdbaselineaf{NO}=[];AUCaf{NO}=[];
    zAVGbaselineaf{1}=mean(mouseztracesmeanpercellaf{i}(1*FR+1:3*FR,:));%Change the number if you want to define the baseline in a different time frame
    zAVGResponseaf{1}=mean(mouseztracesmeanpercellaf{i}(BS*FR+1:(BS+2)*FR,:));
    zStdbaselineaf{1}=std(mouseztracesmeanpercellaf{i}(1*FR+1:3*FR,:));

%%
RCallaf{NO}=[];NCallaf{NO}=[];ICallaf{NO}=[];
     s=size(mouseztracesmeanpercellaf{i},2);
     RC=[];NC=[];IC=[];
     for j=1:s
         if zAVGResponseaf{1}(j)>3*zStdbaselineaf{1}(j)
             RCstd=j;
             RC=cat(1,RC,RCstd);
             clear RCstd
         elseif zAVGResponseaf{1}(j)<-3*zStdbaselineaf{1}(j)
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
    
     RCallaf{1}=RC;
     NCallaf{1}=NC;
     ICallaf{1}=IC;

%% Subtract Mineral oil (MO) responsive cells 
submouseztracesmeanpercellaf=mouseztracesmeanpercellaf;
for i=1:5
    submouseztracesmeanpercellaf{i}(:,RCallaf{1})=mouseztracesmeanpercellaf{i}(:,RCallaf{1})-mouseztracesmeanpercellaf{1}(:,RCallaf{1});    
end

for i=1:5
   submouseztracesmeanpercellaf{i}(:,ICallaf{1}(j,1))=mouseztracesmeanpercellaf{i}(:,ICallaf{10}(j,1))-mouseztracesmeanpercellaf{10}(:,ICallaf{10}(j,1));   
end

%% plot the traces of all the animals before condition
for i=1:size(submouseztracesmeanpercell,2)
figure(i) 
plot(Time4plot,mean(submouseztracesmeanpercell{i},2),'LineWidth',2,'color','k');
hold on
SE = std(submouseztracesmeanpercell{i}')/sqrt(size(submouseztracesmeanpercell{i},2)); %I am not sure about what the SE will be??(:,RCall{i})
    CIplus = mean(submouseztracesmeanpercell{i},2)+(1.96*SE');%(:,RCall{i})
    CIminus = mean(submouseztracesmeanpercell{i},2)-(1.96*SE');%(:,RCall{i})
  shade(Time4plot,CIplus,Time4plot,CIminus,'FillType',[1 2;2 1],'LineStyle','none','FillColor',[0.8 0.8 0.8])%red shadow [0.6,0,0]
    clear CIplus CIminus SE  
xline(0,'LineWidth',1) %on line
%ax.XAxis.FontSize=30
xline(2,'LineWidth',1) % off line % 1 sec odor presentation
ylim([-0.5 1])
%clear mean
end
%% plot the traces of all the animals after condition
for i=1:size(submouseztracesmeanpercellaf,2)
figure(i)    
   
plot(Time4plot,mean(submouseztracesmeanpercellaf{i},2),'LineWidth',2,'color','r');
hold on
SE = std(submouseztracesmeanpercellaf{i}(:,142:end)')/sqrt(size(submouseztracesmeanpercellaf{i}(:,142:end),2)); %I am not sure about what the SE will be??(:,RCall{i})
    CIplus = mean(submouseztracesmeanpercellaf{i},2)+(1.96*SE');%(:,RCall{i})
    CIminus = mean(submouseztracesmeanpercellaf{i}(:,142:end),2)-(1.96*SE');%(:,RCall{i})
  shade(Time4plot,CIplus,Time4plot,CIminus,'FillType',[1 2;2 1],'LineStyle','none','FillColor',[0.6,0,0])%red shadow [0.6,0,0]
    clear CIplus CIminus SE  
xline(0,'LineWidth',1) %on line
xline(2,'LineWidth',1) % off line % 1 sec odor presentation
%ax.XAxis.FontSize=30
ylim([-0.5 1])
%clear mean
end

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
%% Calculate the response of neurons that are activated/inhibited by both Pentanol and Hexanol after condition
Common_RC_af=intersect(RCallaf{2},RCallaf{3});
Common_IC_af=intersect(ICallaf{2},ICallaf{3});
zAVG_CommonRC=cat(1,zAVGResponseaf{2}(:,Common_RC_af), zAVGResponseaf{3}(:,Common_RC_af));
zAVG_CommonIC=cat(1,zAVGResponseaf{2}(:,Common_IC_af), zAVGResponseaf{3}(:,Common_IC_af));
%% Make the traces ready for PCA
PCAmousetracesmeanpercellbf=[];PCAmousetracesmeanpercellaf=[];
for i=2:5
    PCAmousetracesmeanpercellbf=cat(1,PCAmousetracesmeanpercellbf,submouseztracesmeanpercell{i}(:,104:end));
    PCAmousetracesmeanpercellaf=cat(1,PCAmousetracesmeanpercellaf,submouseztracesmeanpercellaf{i}(:,104:end));

end
%%
PCAbfaf_pent_Hex_Oct=[];
PCAbfaf_pent_Hex_Oct=cat(1,PCAbfaf_pent_Hex_Oct,PCAmousetracesmeanpercellbf(1:392,:),PCAmousetracesmeanpercellbf(589:784,:),PCAmousetracesmeanpercellaf(1:392,:),PCAmousetracesmeanpercellaf(589:784,:));
%% Before and after condition (Pentanol vs Hexanol vs Octanol)
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

%% Plot Pentanol Octanol before after
figure();hold on
for ii = 1:6
    if ii==2
        continue
    end
    if ii==5 
       continue
    end
 
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
%% Plot Hexanol Octanol before after
figure();hold on
for ii = 2:6
   
    if ii==4 
       continue
    end
 
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


%% PCA_Before condition four odors
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
   %% PCA_After condition four odors
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
%% Calculate the Euclidean distance before and after condition for each odor
samplingRate = 1;  % Frames per second, adjust based on your data
windowDuration = 1;  % Duration in seconds for the odor presentation window
windowSize = windowDuration * samplingRate;  % Calculate the number of frames in the time window
numTimePoints = 196;  % Total number of time points in the data (adjust based on your data)
numComparisons = numTimePoints - windowSize;  % Number of consecutive time window comparisons

% Initialize arrays to store angles and distances for each neuron, trial, and comparison
cDBetweenSets=zeros(196,8,5);
%anglesBetweenSetsaf = zeros(196,8,5);
distanceBetweenSetsaf = zeros(196,8,5);
%medallaf=zeros(5,900);
% Loop through each trial to calculate the angle and distance
for i = 1:5
    for mouse=1:8
        % Loop through consecutive time windows for each trial
        for t = 1:196  % Loop from the first to the penultimate window
            % Extract calcium activity for the current and previous time windows
            currentWindow = mousenum{1,mouse}.ztemptraceall{i}(t, :,:)-mean(mousenum{1,mouse}.ztemptraceall{i}(1:45, :,:),1);
            previousWindow =mousenum{1,mouse}.ztemptraceallaf{i}(t, :,:)-mean(mousenum{1,mouse}.ztemptraceallaf{i}(1:45, :,:),1);
           
            
            % % Calculate the mean activity vector over the time window
            meanCurrentActivity = mean(currentWindow, 3);  % Mean across time points
            meanPreviousActivity = mean(previousWindow, 3);  
            
            % Calculate the dot product between the two mean activity vectors
            dotProduct = dot(meanPreviousActivity,meanCurrentActivity);
            
            % ------ Euclidean Distance ------
            euclideanDistance = norm(meanCurrentActivity- meanPreviousActivity);  % Euclidean distance
            
            % Calculate the norms (magnitudes) of both activity vectors
            normCurrentActivity = norm(meanCurrentActivity);  % Magnitude of the current activity vector
            normPreviousActivity = norm(meanPreviousActivity);  % Magnitude of the previous activity vector
            
            % Calculate the angle in radians between the two vectors using acos
            %angle = acos(dotProduct / (normCurrentActivity * normPreviousActivity));  % Angle in radians
            cosineDistance = 1-(dotProduct / (normPreviousActivity * normCurrentActivity ));
            
            % Store the angle and Euclidean distance for this comparison
            cDBetweenSets(t,mouse,i)=cosineDistance;
            %anglesBetweenSetsaf(t,mouse,i) = angle;   % Convert to degrees
            distanceBetweenSetsaf(t,mouse,i) = euclideanDistance;
         end
      
    end
end
%% plot the Euclidean distance before and after condition along time for each odor
%for i=1:size(submouseztracesmeanpercell,2)
clrshade=[0.8 0.8 0.8;0,1,0;0.92,0.61,0.61;0.93,0.69,0.13;0.06,1,1];
clr=[0,0,0;0,0.3,0.2;0.64,0.08,0.18;1,1,0;0,0,1];


for i=2:5 
    figure(1) 
    if i==4
        continue
    end

plot(Time4plot,mean(distanceBetweenSetsaf(:,:,i),2),'LineWidth',2,'color',clr(i,:));
hold on
SE = std(distanceBetweenSetsaf(:,:,i)')/sqrt(size(distanceBetweenSetsaf(:,:,i),2)); %I am not sure about what the SE will be??(:,RCall{i})
    CIplus = mean(distanceBetweenSetsaf(:,:,i),2)+(1.96*SE');%(:,RCall{i})
    CIminus = mean(distanceBetweenSetsaf(:,:,i),2)-(1.96*SE');%(:,RCall{i})
  shade(Time4plot,CIplus,Time4plot,CIminus,'FillType',[1 2;2 1],'LineStyle','none','FillColor',clrshade(i,:))%red shadow [0.6,0,0]
    clear CIplus CIminus SE  
xline(0,'LineWidth',1) %on line
%ax.XAxis.FontSize=30
xline(2,'LineWidth',1) % off line % 1 sec odor presentation
% ylim([-0.5 1])
hold on

end
%% Calculate the Euclidean distance from the control odor octanol before condition
% Initialize arrays to store angles and distances for each neuron, trial, and comparison
cDtocontrol=zeros(196,8,4);
%anglesBetweenSetsaf = zeros(196,8,5);
distancetocontrol = zeros(196,8,4);
%medallaf=zeros(5,900);
% Loop through each trial to calculate the angle and distance
for i = 1:4
    for mouse=1:8
        % Loop through consecutive time windows for each trial
        for t = 1:196  % Loop from the first to the penultimate window
            % Extract calcium activity for the current and previous time windows
            currentWindow = mousenum{1,mouse}.ztemptraceall{i}(t, :,:);%-mean(mousenum{1,mouse}.ztemptraceall{i}(1:45, :,:),1);
            previousWindow =mousenum{1,mouse}.ztemptraceall{5}(t, :,:);%-mean(mousenum{1,mouse}.ztemptraceallaf{i}(1:45, :,:),1);
           
            
            % % Calculate the mean activity vector over the time window
            meanCurrentActivity = mean(currentWindow, 3);  % Mean across time points
            meanPreviousActivity = mean(previousWindow, 3);  
            
            % Calculate the dot product between the two mean activity vectors
            dotProduct = dot(meanPreviousActivity,meanCurrentActivity);
            
            % ------ Euclidean Distance ------
            euclideanDistance = norm(meanCurrentActivity- meanPreviousActivity);  % Euclidean distance
            
            % Calculate the norms (magnitudes) of both activity vectors
            normCurrentActivity = norm(meanCurrentActivity);  % Magnitude of the current activity vector
            normPreviousActivity = norm(meanPreviousActivity);  % Magnitude of the previous activity vector
            
            % Calculate the angle in radians between the two vectors using acos
            %angle = acos(dotProduct / (normCurrentActivity * normPreviousActivity));  % Angle in radians
            cosineDistance = 1-(dotProduct / (normPreviousActivity * normCurrentActivity ));
            
            % Store the angle and Euclidean distance for this comparison
            cDtocontrol(t,mouse,i)=cosineDistance;
            %anglesBetweenSetsaf(t,mouse,i) = angle;   % Convert to degrees
            distancetocontrol(t,mouse,i) = euclideanDistance;
         end
      
    end
end
%% plot the Euclidean distance from the control odor octanol before condition along time 
%for i=1:size(submouseztracesmeanpercell,2)
clrshade=[0.8 0.8 0.8;0,1,0;0.92,0.61,0.61;0.93,0.69,0.13;0.06,1,1];
clr=[0,0,0;0,0.3,0.2;0.64,0.08,0.18;1,1,0;0,0,1];


for i=2:4 
    figure(2) 
  

plot(Time4plot,mean(distancetocontrol(:,:,i),2),'LineWidth',2,'color',clr(i,:));
hold on
SE = std(distancetocontrol(:,:,i)')/sqrt(size(distancetocontrol(:,:,i),2)); %I am not sure about what the SE will be??(:,RCall{i})
    CIplus = mean(distancetocontrol(:,:,i),2)+(1.96*SE');%(:,RCall{i})
    CIminus = mean(distancetocontrol(:,:,i),2)-(1.96*SE');%(:,RCall{i})
  shade(Time4plot,CIplus,Time4plot,CIminus,'FillType',[1 2;2 1],'LineStyle','none','FillColor',clrshade(i,:))%red shadow [0.6,0,0]
    clear CIplus CIminus SE  
xline(0,'LineWidth',1) %on line
%ax.XAxis.FontSize=30
xline(2,'LineWidth',1) % off line % 1 sec odor presentation
% ylim([-0.5 1])
hold on

end
%% Calculate the Euclidean distance from the control odor octanol after condition
% Initialize arrays to store angles and distances for each neuron, trial, and comparison
cDtocontrolaf=zeros(196,8,4);
%anglesBetweenSetsaf = zeros(196,8,5);
distancetocontrolaf = zeros(196,8,4);
%medallaf=zeros(5,900);
% Loop through each trial to calculate the angle and distance
for i = 1:4
    for mouse=1:8
        % Loop through consecutive time windows for each trial
        for t = 1:196  % Loop from the first to the penultimate window
            % Extract calcium activity for the current and previous time windows
            currentWindow = mousenum{1,mouse}.ztemptraceallaf{i}(t, :,:);%-mean(mousenum{1,mouse}.ztemptraceall{i}(1:45, :,:),1);
            previousWindow =mousenum{1,mouse}.ztemptraceallaf{5}(t, :,:);%-mean(mousenum{1,mouse}.ztemptraceallaf{i}(1:45, :,:),1);
           
            
            % % Calculate the mean activity vector over the time window
            meanCurrentActivity = mean(currentWindow, 3);  % Mean across time points
            meanPreviousActivity = mean(previousWindow, 3);  
            
            % Calculate the dot product between the two mean activity vectors
            dotProduct = dot(meanPreviousActivity,meanCurrentActivity);
            
            % ------ Euclidean Distance ------
            euclideanDistance = norm(meanCurrentActivity- meanPreviousActivity);  % Euclidean distance
            
            % Calculate the norms (magnitudes) of both activity vectors
            normCurrentActivity = norm(meanCurrentActivity);  % Magnitude of the current activity vector
            normPreviousActivity = norm(meanPreviousActivity);  % Magnitude of the previous activity vector
            
            % Calculate the angle in radians between the two vectors using acos
            %angle = acos(dotProduct / (normCurrentActivity * normPreviousActivity));  % Angle in radians
            cosineDistance = 1-(dotProduct / (normPreviousActivity * normCurrentActivity ));
            
            % Store the angle and Euclidean distance for this comparison
            cDtocontrolaf(t,mouse,i)=cosineDistance;
            %anglesBetweenSetsaf(t,mouse,i) = angle;   % Convert to degrees
            distancetocontrolaf(t,mouse,i) = euclideanDistance;
         end
      
    end
end
%% plot the Euclidean distance from the control odor octanol after condition along time 
%for i=1:size(submouseztracesmeanpercell,2)
%clrshade=[0.8 0.8 0.8;0,1,0;0.92,0.61,0.61;0.93,0.69,0.13;0.06,1,1];
%clr=[0,0,0;0,0.3,0.2;0.64,0.08,0.18;1,1,0;0,0,1];
clr = [0, 0.5, 0;   % Green
       1, 0, 0;   % Red
       1, 0.5, 0]; % Orange
lightenFactor = 0.5;
clr_lighter = clr + (1 - clr) * 0.3;
clrshade= clr + (1 - clr) * 0.5;
clr_shade_lighter=clrshade + (1 - clrshade) * 0.5;

for i=2:4 
    figure(i) 
  

plot(Time4plot,mean(distancetocontrol(:,:,i),2),'LineWidth',2,'color',clr_lighter((i-1),:));
hold on
SE = std(distancetocontrol(:,:,i)')/sqrt(size(distancetocontrol(:,:,i),2)); %I am not sure about what the SE will be??(:,RCall{i})
    CIplus = mean(distancetocontrol(:,:,i),2)+(1.96*SE');%(:,RCall{i})
    CIminus = mean(distancetocontrol(:,:,i),2)-(1.96*SE');%(:,RCall{i})
  shade(Time4plot,CIplus,Time4plot,CIminus,'FillType',[1 2;2 1],'LineStyle','none','FillColor',clr_shade_lighter((i-1),:))%red shadow [0.6,0,0]
    clear CIplus CIminus SE  
xline(0,'LineWidth',1) %on line
%ax.XAxis.FontSize=30
xline(2,'LineWidth',1) % off line % 1 sec odor presentation
ylim([-0.5 6])
hold on
plot(Time4plot,mean(distancetocontrolaf(:,:,i),2),'LineWidth',2,'color',clr((i-1),:));
hold on
SE = std(distancetocontrolaf(:,:,i)')/sqrt(size(distancetocontrol(:,:,i),2)); %I am not sure about what the SE will be??(:,RCall{i})
    CIplus = mean(distancetocontrolaf(:,:,i),2)+(1.96*SE');%(:,RCall{i})
    CIminus = mean(distancetocontrolaf(:,:,i),2)-(1.96*SE');%(:,RCall{i})
  shade(Time4plot,CIplus,Time4plot,CIminus,'FillType',[1 2;2 1],'LineStyle','none','FillColor',clrshade((i-1),:))%red shadow [0.6,0,0]
    clear CIplus CIminus SE  
xline(0,'LineWidth',1) %on line
%ax.XAxis.FontSize=30
xline(2,'LineWidth',1) % off line % 1 sec odor presentation
ylim([-0.5 6])
hold on
end




%% Calculate the Reliability before condition

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

%% Calculate the Reliability after condition

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

%% Save any figure files
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