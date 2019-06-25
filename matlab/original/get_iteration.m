function q = get_iteration(T, ind, deg, plot)
    close all
    q = zeros(deg, 1);
    for i = 1:deg
        vals = SGorthoPolyspk(T, i);
        q(i) = vals(ind);
    end
    
    if plot
        figure
        scatter(q(1:end-1), q(2:end))
        title("Iteration Map for point " + num2str(ind))
        
    end
end