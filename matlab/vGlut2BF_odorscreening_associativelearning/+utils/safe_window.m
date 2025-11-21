% +utils/safe_window.m
function [i0,i1] = safe_window(center, pre, post, N)
i0 = max(1, center-pre);
i1 = min(N, center+post);
end
