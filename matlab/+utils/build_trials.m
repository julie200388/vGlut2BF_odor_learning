% +utils/build_trials.m
function [Trials, tvec] = build_trials(Time, AC, onsets, FR, BS, RS)
nCells = size(AC,2);
winPre = BS*FR; winPost = RS*FR;
len = winPre + winPost + 1;
Trials = NaN(len, nCells, numel(onsets));
for i = 1:numel(onsets)
    [~, idx] = min(abs(Time - onsets(i)));
    [i0,i1] = utils.safe_window(idx, winPre, winPost, size(AC,1));
    seg = AC(i0:i1, :);
    % pad if needed
    if size(seg,1) < len
        padTop = winPre - (idx - i0);
        padBot = (i1 - idx) - winPost;
        Trials(:,:,i) = padarray(seg,[padTop,0],'pre');
        Trials(:,:,i) = padarray(Trials(:,:,i),[padBot,0],'post');
    else
        Trials(:,:,i) = seg;
    end
end
tvec = ((-winPre):winPost)/FR;
end
