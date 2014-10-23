function save_fig(param, label)
    n = 0;
    scope = evalin('base', param);
    save(param, 'scope');
    n = size(scope.signals, 2);
    %n = 1;
    for i = 1:n
        f = figure;
        set(f,'Visible','off');
        set(f, 'Position', [0 0 1000 400]);
        set(f, 'color', 'white');
        p = plot(scope.time, scope.signals(i).values);
        %p = plot(scope);
        
        xlabel('time, $s$', 'Interpreter', 'LaTex');
        %xlabel('Angle error, degree', 'Interpreter', 'LaTex');
        ylabel(label, 'Interpreter', 'LaTex');
        eval(['export_fig -eps -png -m4 ' , param, num2str(i)]);
        close
    end
end