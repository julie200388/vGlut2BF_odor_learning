
clear all

n = input('number of mice: ');
for i = 1:n
    tracevar = uigetfile('*.*');
    mousenum{i} = load(tracevar);
    name1 = strsplit(tracevar, '_');
    name1  = char(name1{1,1});
    mouseID{i} = name1(1:5);
end
%%
dr_mice=zeros(400,5);
for mouse=1:n
    
    kall=[];
for i=1:size(mousenum{1,mouse}.data,2)
k=convertCharsToStrings(mousenum{1,mouse}.data(i).subject);
kall=cat(1,kall,k);
end
s=find(kall=="odor 84_30");

if size(s,1)>1
    drall=[];
    for i=1:size(s,1)
        dr=mousenum{1,mouse}.data(s(i)).dprimeRollingWindow;
        drall=cat(1,drall,dr);
    end
else 
    size(s,1)==1;
drall=mousenum{1,mouse}.data(s(1)).dprimeRollingWindow;
end


dr_mice(1:size(drall,1),mouse)=drall;

end
%%
dr_mice(find(dr_mice==0))=NaN;
