function dom = dominates(a, b)
    % Check if a dominates b
    dom = all(a <= b) && any(a < b);
end