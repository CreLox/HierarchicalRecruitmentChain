% [Stats, PhoR] = LvStats(RepeatNum, SimulationFcn, varargin)
function [Stats, PhoR] = LvStats(RepeatNum, SimulationFcn, varargin)
    tic;
    [~, X] = SimulationFcn(varargin{:});
    TotLv = size(X, 2) - 1;
    Stats = zeros(RepeatNum, TotLv);
    PhoR = zeros(RepeatNum, TotLv);
    for i = 1 : RepeatNum
        [T, X] = SimulationFcn(varargin{:});
        PhoR(i, :) = ActualPhoR(T, X);
        for j = 1 : size(X, 1) - 1
            Stats(i, X(j, end)) = Stats(i, X(j, end)) + T(j + 1) - T(j);
        end
    end
    TotT = sum(Stats, 2);
    for i = 1 : TotLv
        Stats(:, i) = Stats(:, i) ./ TotT;
    end
    toc;
end
