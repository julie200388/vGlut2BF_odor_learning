% +utils/pca_and_plots.m

function pca_and_plots(MeanPerCell, figdir, tag, isAfter, paired)
% stacks odors vertically: [time x cells] for each odor
if nargin < 4, isAfter = false; end     % default: not "after"
if nargin < 5, paired  = false; end     % default: not paired bf+af
% stacks odors vertically: [time x cells] for each odor
X = [];
lbl = [];
for k = 1:numel(MeanPerCell)
    M = MeanPerCell{k};
    if isempty(M), continue; end
    X = [X; M]; %#ok<AGROW>
    lbl = [lbl; k*ones(size(M,1),1)]; %#ok<AGROW>
end
if isempty(X), return; end
CV = cov(X);
[ev, el] = eig(CV);
var = diag(el)/sum(diag(el));
data3 = X*ev(:,end-2:end);
figure('Color','w'); hold on

K = max(lbl);
if paired

cmap = lines(K/2);
else
cmap = lines(K); % base colors for 'before'
end
% Create darker versions for 'after' odors
darkFactor = 0.6; % smaller = darker
cmap_after = cmap * darkFactor;
cmapcombined = [cmap; cmap_after];

% If this is the "after" condition, make colors darker
if isAfter

    cmap2 = cmap_after;

elseif paired
    cmap2= cmapcombined;

else
    cmap2 = cmap;

end
% Combine them vertically

for k = 1:K
    seg = data3(lbl==k,:);
    plot3(seg(:,3), seg(:,2), seg(:,1), '-', 'LineWidth',2, 'Color', cmap2(k,:));
    plot3(seg(1,3), seg(1,2), seg(1,1), 'o','MarkerFaceColor',cmap2(k,:),'MarkerEdgeColor','none');
end
grid on; view([210 30]);
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
legend(arrayfun(@(k)sprintf('odor %d',k),1:K,'uni',0),'Location','northeast','Box','off');

saveas(gcf, fullfile(figdir, sprintf('pca_%s.png',tag)));
close
end
