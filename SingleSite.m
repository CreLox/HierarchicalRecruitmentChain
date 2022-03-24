% [T, X] = SingleSite(TotLv, phorate, dephorate, InitialStates)
function [T, X] = SingleSite(TotLv, phorate, dephorate, InitialStates)
    %% Parameter setting & initialization
    Params.TotLv         = TotLv; % int, >= 2
    
    % P -> P(pho) ([ATP] incorporated)
    Params.PhoRates      = [ones(1, Params.TotLv - 1), 0] * phorate;
    
    % P(pho) -> P
    Params.DephoRates    = [ones(1, Params.TotLv - 1), 0] * dephorate;
    
    % ...-P(n-1,pho) + P(n) -> ...-P(n-1,pho)-P(n)
    % (concentraions of P's incorporated)
    onrate               = 1;
        Params.OnRates   = [0, ones(1, Params.TotLv - 1), 0] * onrate;
        % The last 0 of Params.OnRates is a fake lv to simplify calculation
        % of propensities (see PropensityFcn)
    
    % ...-P(n-1)-P(n)-... -> ...-P(n-1) + P(n)-...
    offrate              = 1;
        Params.OffRates  = [0, ones(1, Params.TotLv - 1)] * offrate;
    
    TLim                 = 1e5;
    EventLim             = 1e6;
    X                    = zeros(EventLim, Params.TotLv + 1);
        % Col 1 to TotLv: phophorylation states of lv #1 - #TotLv:
        Phosphorylated   = 1;
        Dephosphorylated = 0;
        % The phophorylation state of #TotLv/any unoccupied lv is always 0
        % End Col: highest lv occupied (>=1, as lv #1 is always occupied)
        X(1, :)          = InitialStates;
    T                    = zeros(EventLim, 1);
        T(1)             = 0;
    EventCount           = 1;
    
    %% Discrete-event simulation
    while true
        if EventCount + 1 > EventLim
            % disp('Number of events exceeds the limit.');
            break;
        end
        
        % Calculate propensities & simulate the next event
        Propensities = PropensityFcn(X(EventCount, :), Params);
        PropensitySum = sum(Propensities);
        r = rand(1, 2);
        DeltaT = - log(r(1)) / PropensitySum;
        EventID = find((cumsum(Propensities) >= r(2) * PropensitySum), 1);
        
        % Update the clock and carry out the event
        if T(EventCount) + DeltaT > TLim
            % disp('Time exceeds the limit.');
            T(EventCount + 1) = TLim;
            X(EventCount + 1, :) = X(EventCount, :);
            EventCount = EventCount + 1;
            break;
        else
            T(EventCount + 1) = T(EventCount) + DeltaT;
            switch EventID
                case num2cell(1 : Params.TotLv),
                    % Phosphorylation
                    X(EventCount + 1, 1 : end) = X(EventCount, 1 : end);
                    X(EventCount + 1, EventID) = Phosphorylated;
                case num2cell(Params.TotLv + 1 : 2 * Params.TotLv),
                    % Dephosphorylation
                    X(EventCount + 1, 1 : end) = X(EventCount, 1 : end);
                    X(EventCount + 1, EventID - Params.TotLv) = Dephosphorylated;
                case num2cell(2 * Params.TotLv + 1 : 3 * Params.TotLv),
                    % Binding
                    X(EventCount + 1, end) = X(EventCount, end) + 1;
                    X(EventCount + 1, 1 : Params.TotLv) = X(EventCount, 1 : Params.TotLv);
                case num2cell(3 * Params.TotLv + 1 : 4 * Params.TotLv),
                    % Falling off
                    X(EventCount + 1, end) = EventID - 3 * Params.TotLv - 1;
                    X(EventCount + 1, 1 : (EventID - 3 * Params.TotLv - 1)) = ...
                        X(EventCount, 1 : (EventID - 3 * Params.TotLv - 1));
            end
            EventCount = EventCount + 1;
        end
    end
    X = X(1 : EventCount, :);
    T = T(1 : EventCount);
end

function Propensities = PropensityFcn(States, Params)
    PhoStates = States(1 : Params.TotLv);
    HighestOccupied = States(end);
    
    PhoProp = [~PhoStates(1 : HighestOccupied), ...
        zeros(1, Params.TotLv - HighestOccupied)] .* Params.PhoRates;
    
    DephoProp = PhoStates .* Params.DephoRates;
    
    OnProp = zeros(1, Params.TotLv + 1); % The last element is a fake lv
    OnProp(HighestOccupied + 1) = PhoStates(HighestOccupied) * ...
        Params.OnRates(HighestOccupied + 1);
    
    OffProp = zeros(1, Params.TotLv);
    OffProp(2 : HighestOccupied) = (~PhoStates(1 : HighestOccupied - 1)) .* ...
        Params.OffRates(2 : HighestOccupied);
    
    Propensities = [PhoProp, DephoProp, OnProp(1 : end - 1), OffProp];
end
