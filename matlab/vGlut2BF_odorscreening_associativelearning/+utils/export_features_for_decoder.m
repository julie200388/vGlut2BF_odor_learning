% +utils/export_features_for_decoder.m
function features = export_features_for_decoder(Z_byOdor, winSec, FR)
% stacks trials rows; last column = odor label
w0 = winSec(1)*FR + 1; w1 = winSec(2)*FR;
features = [];
labels   = [];
for k = 1:numel(Z_byOdor)
    Z = Z_byOdor{k};
    if isempty(Z), continue; end
    % reshape time window x cells x trial -> trials x (cells*frames)
    X = permute(Z(w0:w1,:,:), [3 2 1]); % trials x cells x frames
    X = reshape(X, size(X,1), []);      % trials x (cells*frames)
    features = [features; X]; %#ok<AGROW>
    labels   = [labels; k*ones(size(X,1),1)]; %#ok<AGROW>
end
features = [features labels];
end
