% +utils/zscore_trials.m
function Z = zscore_trials(Trials, FR, BS)
% baseline window: frames FR+1 : BS*FR (exclude the first second if desired)
b0 = FR+1; b1 = BS*FR;
Z = Trials;
for i = 1:size(Trials,3)
    base = Trials(b0:b1,:,i);         % time x cells
    mu   = mean(base,1);              % 1 x cells
    sd   = std(base,0,1);             % 1 x cells
    sd(sd==0) = 1;
    Z(:,:,i) = (Trials(:,:,i) - mu)./sd;
end
Z(isnan(Z)) = 0;
end
