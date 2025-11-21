% +utils/split_sessions_by_time.m
function idx = split_sessions_by_time(Time, gapThresh)
% SPLIT_SESSIONS_BY_TIME
% Find indices where the recording time "wraps" to a new session
% (i.e., large negative jump in Time).
%
% INPUTS
%   Time      : N x 1 vector of time stamps (seconds)
%   gapThresh : scalar; a jump < -gapThresh is considered a session boundary
%
% OUTPUT
%   idx       : vector of indices where a new session starts (the last
%               sample of the previous session is at idx, and the first
%               sample of the next session is at idx+1).
%
% Example:
%   - 2 days  -> one wrap index: [D2start]
%   - 3 days  -> two wrap indices: [D2start, D3start]

if nargin < 2
    gapThresh = 1.0;  % default 1 second, adjust if needed
end

Time = Time(:);          % ensure column
dT   = diff(Time);       % size N-1

% Find *all* wraps where time goes backwards by more than gapThresh
wraps = find(dT < -gapThresh);

if isempty(wraps)
    error('split_sessions_by_time:NoWrapsFound', ...
        'Could not find any session splits (no negative jumps < -%.3f).', gapThresh);
end

idx = wraps;   % return all wrap indices as a vector
end
