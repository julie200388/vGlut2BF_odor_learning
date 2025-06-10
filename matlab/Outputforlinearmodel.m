%% Find MINERAL OIL responsive cells using 3*STD threshold

% Calculate baseline and response statistics
zAVGbaseline = ztracesmeanpercell{1}(1*FR+1:3*FR, :);
zAVGResponse = mean(ztracesmeanpercell{1}(BS*FR+1:(BS+2)*FR, :), 1);
zStdbaseline = std(zAVGbaseline, 1);

zAVGbaselineaf = ztracesmeanpercellaf{1}(1*FR+1:3*FR, :);
zAVGResponseaf = mean(ztracesmeanpercellaf{1}(BS*FR+1:(BS+2)*FR, :), 1);
zStdbaselineaf = std(zAVGbaselineaf, 1);

%% Classify cells for "before" condition
RC = []; IC = []; NC = [];
for j = 1:size(ztracesmeanpercell{1}, 2)
    if zAVGResponse(j) > 3 * zStdbaseline(j)
        RC = [RC; j];
    elseif zAVGResponse(j) < -3 * zStdbaseline(j)
        IC = [IC; j];
    else
        NC = [NC; j];
    end
end

if isempty(RC), RC = 0; end
if isempty(IC), IC = 0; end
if isempty(NC), NC = 0; end

RCbf = RC; ICbf = IC; NCbf = NC;

%% Classify cells for "after" condition
RC = []; IC = []; NC = [];
for j = 1:size(ztracesmeanpercellaf{1}, 2)
    if zAVGResponseaf(j) > 3 * zStdbaselineaf(j)
        RC = [RC; j];
    elseif zAVGResponseaf(j) < -3 * zStdbaselineaf(j)
        IC = [IC; j];
    else
        NC = [NC; j];
    end
end

if isempty(RC), RC = 0; end
if isempty(IC), IC = 0; end
if isempty(NC), NC = 0; end

RCaf = RC; ICaf = IC; NCaf = NC;

%% Subtract MO (mineral oil) response from all conditions
for i = 1:5
    if RCbf ~= 0
        ztemptraceall{i}(:, RCbf, :) = ztemptraceall{i}(:, RCbf, :) - ztemptraceall{1}(:, RCbf, :);
    end
    if RCaf ~= 0
        ztemptraceallaf{i}(:, RCaf, :) = ztemptraceallaf{i}(:, RCaf, :) - ztemptraceallaf{1}(:, RCaf, :);
    end
    if ICbf ~= 0
        ztemptraceall{i}(:, ICbf, :) = ztemptraceall{i}(:, ICbf, :) - ztemptraceall{1}(:, ICbf, :);
    end
    if ICaf ~= 0
        ztemptraceallaf{i}(:, ICaf, :) = ztemptraceallaf{i}(:, ICaf, :) - ztemptraceallaf{1}(:, ICaf, :);
    end
end

%% Extract traces 2 seconds after odor onset for linear model input
traces4model2secbf = [];
traces4model2secaf = [];
NO = size(ztemptraceall, 2);

for i = 1:NO
    featureall = [];
    featureallaf = [];
    for j = 1:10
        feature = reshape(ztemptraceall{i}(46:75, :, j), 1, NAC * 45);
        featureaf = reshape(ztemptraceallaf{i}(46:75, :, j), 1, NAC * 45);
        featureall = [featureall; feature];
        featureallaf = [featureallaf; featureaf];
    end
    traces4model2secbf = [traces4model2secbf; featureall];
    traces4model2secaf = [traces4model2secaf; featureallaf];
    odor4feature(10 * (i-1) + 1 : 10 * i, 1) = i;
end

traces4model2secbf = [traces4model2secbf odor4feature];
traces4model2secaf = [traces4model2secaf odor4feature];

% Save if needed
% writematrix(traces4model2secbf,'072624_A134R1_tracesfeature_odorbfcondition_2SECAFODORONSET_v3.csv')
% writematrix(traces4model2secaf,'080324_A134R1_tracesfeature_odorafcondition_2SECAFODORONSET_v3.csv')