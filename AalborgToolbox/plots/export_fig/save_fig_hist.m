function save_fig_hist(param, label)
    n = 0;
    scope = evalin('base', param);
    %save(param, 'scope');
    n = size(scope.signals, 2);
    %n = 1;
    for i = 1:n
        f = figure;
        set(f,'Visible','off');
        set(f, 'Position', [0 0 1000 400]);
        set(f, 'color', 'white');
        hist = histc(scope.signals(i).values, 1:1:180)
        p = plot(hist);
        %p = plot(scope);
        
        %xlabel('time, $s$', 'Interpreter', 'LaTex');
        xlabel(label, 'Interpreter', 'LaTex');
        ylabel('histogram counts, #', 'Interpreter', 'LaTex');
        eval(['export_fig -eps -png -m4 ' , param, '_hist', num2str(i)]);
        close
    end
end