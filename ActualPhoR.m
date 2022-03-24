% r = ActualPhoR(T, X)
function r = ActualPhoR(T, X)
    TotLv = size(X, 2) - 1;
    r = NaN(1, TotLv);
    for i = 1 : TotLv
        RowIdx = find(X(1 : end - 1, end)' >= i);
        if ~isempty(RowIdx)
            TotPhoTime = 0;
            TotDephoTime = 0;
            for j = RowIdx
                if (X(j, i) == 1)
                    TotPhoTime = TotPhoTime + T(j + 1) - T(j);
                else
                    TotDephoTime = TotDephoTime + T(j + 1) - T(j);
                end            
            end
            r(i) = TotPhoTime / (TotPhoTime + TotDephoTime);
        end
    end
end
