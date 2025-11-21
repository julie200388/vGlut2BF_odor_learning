% +utils/load_dio_and_onsets.m
function onsets = load_dio_and_onsets(csvPath, thr)
T = readtable(csvPath);
if size(T,2) >= 2
    T(:,2) = [];  % mirror of your original
end
if height(T) >= 32
    T(1:32,:) = [];
end
A = table2array(T);
bin = double(A(:,2) > thr);
d   = diff(bin);
onsets = A(find(d==1)+1, 1);
end
