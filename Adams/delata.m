function d = delata(f, j, n_idx)
% Help function to caclulate deferences fn and fn-1
    if j == 0
        d = f(n_idx);
    else
        d = delata(f, j - 1, n_idx) - delata(f, j - 1, n_idx - 1);
    end
end