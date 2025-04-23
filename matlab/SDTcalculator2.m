function [SDTout, Dout, f2] = SDTcalculator(ConcatenatedSessions)


together = ConcatenatedSessions;
    
    
    %%% crop to completed trials
    idxend = find(together(:,1) == 0);
    together(idxend,:) = [];
    idxnan = find(isnan(together(:,1))==1);
    together(idxnan,:) = [];
 
    if length(together) > 49
    
    %%% For calculating SDT measures
    %%% d' = sensitivity
    %%% criterion = bias
    
    Gtrials = length(find(together(:,2) == 1))+1;
    NGtrials = length(find(together(:,2) == 0))+1;
    FAtrials = length(find(together(:,3)==-1))+.5;
    Htrials = length(find(together(:,3)==1))+.5;
    
    FA = (FAtrials/NGtrials);
    HR = (Htrials/Gtrials);
    
    dprime = icdf('normal',(HR),0,1)-icdf('normal',(FA),0,1)
    criterion = -1*(icdf('normal',HR,0,1)+icdf('normal',FA,0,1))/2
    
    
    
    %dlmwrite('SDTvalues.txt',[length(together) dprime criterion],'delimiter','\t','-append')
    
    %%% For plotting accuracy by block
    % y = Nvec';
    % x = [20:20:(length(Nvec)*20)];
    % success = figure(1);
    % plot(x,y,'-ro','MarkerSize',10,...
    %     'MarkerEdgeColor','red',...
    %     'MarkerFaceColor',[1 .6 .6])
    % axis([0 max(x) 0 1])
    % yline(.85,'LineStyle','--','Color','g','LineWidth',1)
    % yline(.5,'LineStyle',':','LineWidth',1)
    % x0=200;
    % y0=200;
    % width=600;
    % height=200;
    % set(gcf,'position',[x0,y0,width,height])
    
    %%% Calculating rolling window d'
    window = 40;
    stopidx = length(together)-window;
    
    dprimeRW = [];
    for i = 1:stopidx
        
        GtrialsRW = length(find(together(i:i+window,2) == 1))+1;
        NGtrialsRW = length(find(together(i:i+window,2) == 0))+1;
        FAtrialsRW = length(find(together(i:i+window,3)==-1))+.5;
        HtrialsRW = length(find(together(i:i+window,3)==1))+.5;
        FARW = (FAtrialsRW/NGtrialsRW);
        HRRW = (HtrialsRW/GtrialsRW);
        
        dprimeRW(i) = icdf('normal',HRRW,0,1)-icdf('normal',FARW,0,1);
        
        clear GtrialsRW NGtrialsRW FAtrialsRW HtrialsRW FARW HRRW
    end
    clear i
    
    %%% For plotting rolling window d'
    
    DPmax = icdf('normal',((window-1)/window),0,1)-icdf('normal',(1/window),0,1);
    DPmax = DPmax+(DPmax*.2);
    DPmin = -1*DPmax;
    f2 = figure(2);
    %if dprimeRW > 0
    plot(dprimeRW,'-ro','MarkerSize',5,...
        'MarkerEdgeColor','red',...
        'MarkerFaceColor',[1 .6 .6])
    axis([0 length(dprimeRW) DPmin DPmax])
    yline(2.0729,'LineStyle','--','Color','g','LineWidth',1)
    %yline(0,'LineStyle',':','LineWidth',1)
    x0=200;
    y0=200;
    width=600;
    height=200;
    set(gcf,'position',[x0,y0,width,height])
    %end
    
    
    %dfn = sprintf('%s%s',medpcfn,'_dpdt.txt');
    %dlmwrite(dfn,dprimeRW');
    
    %%% Outputs
    %fn = string(filename(end-30:end-3));
    
    SDTout = [length(together), dprime, criterion];
    Dout = dprimeRW;
    
else
    
    f2 = figure(1)
    SDTout  = [NaN,NaN,NaN,NaN];
    Dout = [];
end
end


