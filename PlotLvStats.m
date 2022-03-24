% PlotLvStats(Stats, phorate, dephorate)
function PlotLvStats(Stats, phorate, dephorate)
    TotLv = size(Stats, 2);
    Predicted = zeros(TotLv, 1);
        for i = 1 : (TotLv - 1)
            Predicted(i) = i * dephorate * phorate ^ (i - 1);
            for j = 1 : i
                Predicted(i) = Predicted(i) / (j * dephorate + phorate);
            end
        end
        Predicted(TotLv) = phorate ^ (TotLv - 1);
        for j = 1 : (TotLv - 1)
            Predicted(TotLv) = Predicted(TotLv) / (j * dephorate + phorate);
        end
    hold on;
    bar_handle = bar(1 : TotLv, [mean(Stats)', Predicted], 'grouped', ...
        'BarWidth', 1, 'LineWidth', 2);
        bar_handle(1).FaceColor = 'w';
        bar_handle(2).FaceColor = 'k';
    errorbar((1 : TotLv) - 0.15, mean(Stats), std(Stats), '.k', 'LineWidth', 2);
    set(gca, 'XTick', 1 : TotLv, 'FontSize', 14, 'LineWidth', 2);
    xlabel('Higest Lv', 'FontSize', 16);
    ylabel('Time frac.', 'FontSize', 16);
    ylim([0, 1]);
    xlim([0, 8]);
    hold off;
end
