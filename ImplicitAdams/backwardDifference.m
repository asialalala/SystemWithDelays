function d = backwardDifference(f, j, n_idx)
% Function to calculate backward difference
    if j == 0
        d = f(n_idx);
    else
        d = backwardDifference(f, j - 1, n_idx) - backwardDifference(f, j - 1, n_idx - 1);
    end
end