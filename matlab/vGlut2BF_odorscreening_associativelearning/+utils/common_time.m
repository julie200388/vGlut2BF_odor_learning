function t = common_time(A, which)
% Return minimal-length time vector across animals for bf/af
len = inf; t_candidate = [];
for a = 1:numel(A)
    if strcmp(which,'bf'), tv = A(a).tvec_bf(:); else, tv = A(a).tvec_af(:); end
    if numel(tv) < len
        len = numel(tv);
        t_candidate = tv;
    end
end
t = t_candidate(1:len);
end
