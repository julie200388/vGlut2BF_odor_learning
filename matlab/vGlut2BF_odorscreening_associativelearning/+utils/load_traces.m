% +utils/load_traces.m
function tbl = load_traces(p)
assert(isfile(p), 'Trace file not found: %s', p);
tbl = readtable(p, 'TextType','string');
end
