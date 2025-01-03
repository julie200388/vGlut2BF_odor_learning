%% gonogoaccuracy
clear all

[medpcfn,medpcpth] = uigetfile('2022-05-21_11h06m.Subject A40 Nt');
medpc =  sprintf('%s%s',medpcpth,medpcfn);
[medpcfile] = importmedpc(medpc);
% 
% 
for i = 1:length(medpcfile.textdata)
medpctext(i) = string(medpcfile.textdata(i));
end

medpctext = medpctext';
Iidx = find(medpctext == 'I:');
Jidx = find(medpctext == 'J:');
Uidx = find(medpctext == 'U:');
Kidx = find(medpctext == 'K:');
Midx = find(medpctext == 'M:');
Nidx = find(medpctext == 'N:');
rowdif = Jidx-Iidx-1;

Istart = 1;
Jstart = Istart+rowdif+1;
Kstart = Kidx-Iidx;
Mstart = Kstart+(ceil(rowdif/20))+rowdif+1;
Nstart = Mstart+rowdif+1;
Ustart = length(medpcfile.data)-rowdif;

Imat = medpcfile.data(Istart:Istart+rowdif-1,:); %Odor presented 1=go
Jmat = medpcfile.data(Jstart:Jstart+rowdif-1,:); % Duration in odor port
Kmat = medpcfile.data(Kstart+1:Kstart+rowdif,:); %
Umat = medpcfile.data(Ustart+1:Ustart+rowdif,:); %Trial initiation time
MMat = medpcfile.data(Mstart+1:Mstart+rowdif,:); %Reward port entry time in hits and fa
Nmat = medpcfile.data(Nstart+1:Nstart+(ceil(rowdif/20)),:); %Accuracy by block
% Umatsec = Umat/100+(on(1)-Umat(1)/100);
% Mmatsec = MMat/100+(on(1)-Umat(1)/100);
% Umatseczero=find(Umatsec==on(1)-Umat(1)/100);
% Umatsec(Umatseczero)=0;
% Mmatseczero=find(Mmatsec==on(1)-Umat(1)/100);
% Mmatsec(Mmatseczero)=0;

RD = MMat-Umat;
RD(RD < 0) = NaN;
RD(RD == 0) = NaN;

rewarddelay = nanmean(nanmean(RD));
Ivec = reshape(Imat',(size(Imat,1)*size(Imat,2)),1);
Kvec = reshape(Kmat',(size(Kmat,1)*size(Kmat,2)),1);
%Uvec = reshape(Umatsec',(size(Umat,1)*size(Umat,2)),1);
Nvec = reshape(Nmat',(size(Nmat,1)*size(Nmat,2)),1);
Nvec(isnan(Nvec)) = [];
Nvec(Nvec == 0) = [];
%Mvec = reshape(Mmatsec',(size(MMat,1)*size(MMat,2)),1);
Jvec = reshape(Jmat',(size(Jmat,1)*size(Jmat,2)),1);



%%% "together" is variable built from MEDPC data where column 1 is trial
%%% times in seconds, column 2 is go or no go, and column 3 is the success
%%% or failure (0 = correct reject or miss, 1 = hit, -1 = false alarm)

%together(:,1) = Uvec; % Trial initiation time (s)
together(:,2) = Ivec; % Odor presentation: 1 = Go, 0 = No Go
together(:,3) = Kvec; % Result: 1 = Hit, 0 = CR or Miss, -1 = false alarm
%together(:,4) = Mvec; % Time of reward port beam break (s)
%together(:,5) = together(:,4) - together(:,1); % Latency from initiation to reward port beam break
together(:,6) = Jvec;
idx(:) = find(together(:,5) < 0);
together(idx,5) = NaN;

x = Nvec';
y = [20:20:(length(Nvec)*20)];
success = figure(4);

plot(x,y,'-ro','MarkerSize',10,...
    'MarkerEdgeColor','red',...
    'MarkerFaceColor',[1 .6 .6])
axis([0 1 0 max(y)])
axis ij
xline(.85,'LineStyle','--','Color','g','LineWidth',1)
xline(.5,'LineStyle',':','LineWidth',1)
x0=20;
y0=20;
width=200;
height=600;
set(gcf,'position',[x0,y0,width,height])