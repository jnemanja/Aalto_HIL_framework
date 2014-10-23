function [output] = Bdot(input)

    K = 15000;

    persistent Bprev;
    if isempty(Bprev)
        Bprev = [0 0 0];
    end
    
    Bx = input(11);
    By = input(12);
    Bz = input(13);

    output = [-K*[Bx-Bprev(1) By-Bprev(2) Bz-Bprev(3)] 0 0 0 0 0 0];
    
    Bprev = [Bx By Bz];
    
end




