%% Analyzing odor screening dataset
% Process the calcium data of individual mouse

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
%% settining up the TTL pulse file for odor presentationtime
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
%% set up the oder of odor presentations
[fn3,pth3] = uigetfile('*.*');
odorsfn = sprintf('%s%s',pth3,fn3);
odors = csvread(odorsfn);
%% get average traces before condition
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
%%Save individual animal matfile in the same folder
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
%% Combine all the animals
mouseztracesmeanpercell{k}=[];
for i=1:k
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
NO=size(mouseztracesmeanpercellaf,2);
zAVGbaseline{NO}=[];zAVGResponse{NO}=[];zStdbaseline{NO}=[];AUCabs{NO}=[]
    
    zAVGbaseline{10}=mean(mouseztracesmeanpercell{i}(1*FR+1:3*FR,:));%Change the number if you want to define the baseline in a different time frame
    zAVGResponse{10}=mean(mouseztracesmeanpercell{i}(BS*FR+1:(BS+2)*FR,:));
   
    AUCabs{10}=sum(abs(mouseztracesmeanpercell{i}(BS*FR+1:(BS+2)*FR,:)));
    zStdbaseline{10}=std(mouseztracesmeanpercell{i}(1*FR+1:3*FR,:));

%%
RCall{NO}=[];NCall{NO}=[];ICall{NO}=[];

     
     RC=[];NC=[];IC=[];
     for j=1:size(mouseztracesmeanpercell{10},2)
         if zAVGResponse{i}(j)>3*zStdbaseline{10}(j)
             RCstd=j;
             RC=cat(1,RC,RCstd);
             clear RCstd
         elseif zAVGResponse{i}(j)<-3*zStdbaseline{10}(j)
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
    
     RCall{10}=RC;
     NCall{10}=NC;
     ICall{10}=IC;


%% Subtract Mineral oil (MO) responsive cells Ch10 is MO
submouseztracesmeanpercell=mouseztracesmeanpercell;submouseztracesmeanpercellaf=mouseztracesmeanpercellaf;
for i=1:k
    for j=1:size(RCall{10},1)
        if RCall{10}(j,1)~=0
           submouseztracesmeanpercell{i}(:,RCall{10}(j,1))=mouseztracesmeanpercell{i}(:,RCall{10}(j,1))-mouseztracesmeanpercell{1}(:,RCall{10}(j,1));
      
        else
        submouseztracesmeanpercell{i}(:,RCall{10}(j,1))=mouseztracesmeanpercell{i}(:,RCall{10}(j,1));
        end
    end
 
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
%% plot the traces of all the animals (Excite, inhibited, non responsive cells)
for i=1:size(submouseztracesmeanpercell,2)
figure(i) 
plot(Time4plot,mean(submouseztracesmeanpercell{i}(:,RCall{i}),2),'LineWidth',2,'color','b');
hold on
SE = std(submouseztracesmeanpercell{i}(:,RCall{i})')/sqrt(size(submouseztracesmeanpercell{i}(:,RCall{i}),2)); %I am not sure about what the SE will be??(:,RCall{i})
    CIplus = mean(submouseztracesmeanpercell{i}(:,RCall{i}),2)+(1.96*SE');%(:,RCall{i})
    CIminus = mean(submouseztracesmeanpercell{i}(:,RCall{i}),2)-(1.96*SE');%(:,RCall{i})
  shade(Time4plot,CIplus,Time4plot,CIminus,'FillType',[1 2;2 1],'LineStyle','none','FillColor',[0.678, 0.847, 0.902])%red shadow [0.6,0,0]
    clear CIplus CIminus SE  
xline(0,'LineWidth',1) %on line
%ax.XAxis.FontSize=30
xline(2,'LineWidth',1) % off line % 1 sec odor presentation
ylim([-0.5 1])
%clear mean

plot(Time4plot,mean(submouseztracesmeanpercell{i}(:,NCall{i}),2),'LineWidth',2,'color','k');
hold on
SE = std(submouseztracesmeanpercell{i}(:,NCall{i})')/sqrt(size(submouseztracesmeanpercell{i}(:,NCall{i}),2)); %I am not sure about what the SE will be??(:,RCall{i})
    CIplus = mean(submouseztracesmeanpercell{i}(:,NCall{i}),2)+(1.96*SE');%(:,RCall{i})
    CIminus = mean(submouseztracesmeanpercell{i}(:,NCall{i}),2)-(1.96*SE');%(:,RCall{i})
  shade(Time4plot,CIplus,Time4plot,CIminus,'FillType',[1 2;2 1],'LineStyle','none','FillColor',[0.8 0.8 0.8])%red shadow [0.6,0,0]
    clear CIplus CIminus SE  
xline(0,'LineWidth',1) %on line
%ax.XAxis.FontSize=30
xline(2,'LineWidth',1) % off line % 1 sec odor presentation
ylim([-0.5 1])
%clear mean

plot(Time4plot,mean(submouseztracesmeanpercell{i}(:,ICall{i}),2),'LineWidth',2,'color','g');
hold on
SE = std(submouseztracesmeanpercell{i}(:,ICall{i})')/sqrt(size(submouseztracesmeanpercell{i}(:,ICall{i}),2)); %I am not sure about what the SE will be??(:,RCall{i})
    CIplus = mean(submouseztracesmeanpercell{i}(:,ICall{i}),2)+(1.96*SE');%(:,RCall{i})
    CIminus = mean(submouseztracesmeanpercell{i}(:,ICall{i}),2)-(1.96*SE');%(:,RCall{i})
  shade(Time4plot,CIplus,Time4plot,CIminus,'FillType',[1 2;2 1],'LineStyle','none','FillColor',[0.564, 0.933, 0.564])%red shadow [0.6,0,0]
    clear CIplus CIminus SE  
xline(0,'LineWidth',1) %on line
%ax.XAxis.FontSize=30
xline(2,'LineWidth',1) % off line % 1 sec odor presentation
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