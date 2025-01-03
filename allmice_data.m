%% plot all the animals
clear all
%close all

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
mouseztracesmeanpercell2{k}=[];
for i=1:k
    for mouse=1:6
        if n==1
       
        
        mouseztracesmeanpercell2{i}=mean(mousenum{1,mouse}.ztemptraceall{i},3);
        end
        
        if n>1
      
        mouseztracesmeanpercell2{i}=cat(2,mouseztracesmeanpercell2{i},mean(mousenum{1,mouse}.ztemptraceall{i},3));  
        end
    end
end
%%
mouseztracesmeanpercellaf{k}=[];
for i=1:k
    for mouse=1:n
        if n==1
        mouseztracesmeanpercellaf{i}=mousenum{1,mouse}.ztracesmeanpercellaf{i};
        end
        
        if n>1
    
        mouseztracesmeanpercellaf{i}=cat(2,mouseztracesmeanpercellaf{i},mousenum{1,mouse}.ztracesmeanpercellaf{i});  
        end
    end
end
%% Combine all the animals_every trials
mouseztemptraces2{k}=[];
for i=1:k
    for mouse=1:n
        if n==1
        mouseztemptraces2{i}=mousenum{1,mouse}.ztemptraceall{i};
        end
        
        if n>1
    
        mouseztemptraces2{i}=cat(2,mouseztemptraces2{i},mousenum{1,mouse}.ztemptraceall{i});  
        end
    end
end
%%
mouseztemptracesaf{k}=[];
for i=1:k
    for mouse=1:n
        if n==1
        mouseztemptracesaf{i}=mousenum{1,mouse}.ztemptraceallaf{i};
        end
        
        if n>1
    
        mouseztemptracesaf{i}=cat(2,mouseztemptracesaf{i},mousenum{1,mouse}.ztemptraceallaf{i});  
        end
    end
end
%% try this for plot the first 5 trials data
mouseztemptraces_5trilas{5}=[];
for i=1:5
    mouseztemptraces_5trilas{i}=cat(2,mouseztemptraces{i},mouseztemptraces{i}(:,:,1:5));
end
%% try this for plot the first 5 trials data
for i=1:5
    mouseztracesmeanpercell{i}=mean(mouseztemptraces_5trilas{i},3);
    mouseztracesmeanpercellaf{i}=mean(mouseztemptracesaf{i}(:,:,1:5),3);
end