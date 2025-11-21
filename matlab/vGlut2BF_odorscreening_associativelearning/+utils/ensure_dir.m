% +utils/ensure_dir.m
function ensure_dir(p)
if ~exist(p,'dir'), mkdir(p); end
end
