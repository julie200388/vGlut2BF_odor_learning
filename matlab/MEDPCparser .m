function [parsedMEDPCoutput] = MEDPCparser()

[medpcfile] = importmedpc('');

for i = 1:size(medpcfile.textdata)
    medpctext(i) = string(medpcfile.textdata(i));
end


%%
% Open the file
fileID = fopen('2022-05-19_07h46m.Subject A38 Nt', 'r');

% Check if file opened successfully
if fileID == -1
    error('File cannot be opened');
end

% Read each line of the file
line = fgetl(fileID);
while ischar(line)
    disp(line)  % or process the line here
    line = fgetl(fileID);  % Read next line
end
%%
% Open the file
fileID = fopen('2022-05-19_07h46m.Subject A38 Nt.txt', 'r');

% Read data according to format (e.g., '%f' for floating-point numbers)
data = fscanf(fileID, '%f');  % Adjust the format specifier as needed

% Close the file
fclose(fileID);

% Display the data
disp(data);
%%
%medpcfn = medpctext(1);

medpctext = medpctext';
Iidx = find(medpctext == 'I:')-3;
Jidx = find(medpctext == 'J:')-3;
Uidx = find(medpctext == 'U:')-3;
Kidx = find(medpctext == 'K:')-3;
Midx = find(medpctext == 'M:')-3;
Nidx = find(medpctext == 'N:')-3;
rowdif = Jidx-Iidx-1;

Istart = 1;
Jstart = Istart+rowdif+1;
Kstart = Kidx-Iidx;
Mstart = Kstart+(ceil(rowdif/20))+rowdif+1;
Nstart = Mstart+rowdif+1;
Ustart = length(medpcfile.data)-rowdif;

Imat = medpcfile.data(Istart:Istart+rowdif-1,:); % Odor presented 1=go
Jmat = medpcfile.data(Jstart:Jstart+rowdif-1,:); % Duration in odor port
Kmat = medpcfile.data(Kstart+1:Kstart+rowdif,:); % reward array
Umat = medpcfile.data(Ustart+1:Ustart+rowdif,:); % Trial initiation time
MMat = medpcfile.data(Mstart+1:Mstart+rowdif,:); %Reward port entry time in hits and fa
Nmat = medpcfile.data(Nstart+1:Nstart+(ceil(rowdif/20)),:); %Accuracy by block
Umatsec = Umat/100;
Mmatsec = MMat/100;

RD = MMat-Umat;
RD(RD < 0) = NaN;
RD(RD == 0) = NaN;

rewarddelay = nanmean(nanmean(RD));

Ivec = reshape(Imat',(size(Imat,1)*size(Imat,2)),1);
Kvec = reshape(Kmat',(size(Kmat,1)*size(Kmat,2)),1);
Uvec = reshape(Umatsec',(size(Umat,1)*size(Umat,2)),1);
Nvec = reshape(Nmat',(size(Nmat,1)*size(Nmat,2)),1);
Nvec(isnan(Nvec)) = [];
Nvec(Nvec == 0) = [];
Mvec = reshape(Mmatsec',(size(MMat,1)*size(MMat,2)),1);
Jvec = reshape(Jmat',(size(Jmat,1)*size(Jmat,2)),1);

%%% "together" is variable built from MEDPC data where column 1 is trial
%%% times in seconds, column 2 is go or no go, and column 3 is the success
%%% or failure (0 = correct reject or miss, 1 = hit, -1 = false alarm)

together(:,1) = Uvec; % Trial initiation time (s)
together(:,2) = Ivec; % Odor presentation: 1 = Go, 0 = No Go
together(:,3) = Kvec; % Result: 1 = Hit, 0 = CR or Miss, -1 = false alarm
together(:,4) = Mvec; % Time of reward port beam break (s)
together(:,5) = together(:,4) - together(:,1); % Latency from initiation to reward port beam break
together(:,6) = Jvec; % Duration in odor port (i.e. crude reaction time)
idx(:) = find(together(:,5) < 0);
together(idx,5) = NaN;

parsedMEDPCoutput = together;


