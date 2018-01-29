function idx_0 = centromere(idx, thr)
d = diff(idx);
A = find(d == 1);
B = find(d == -1);
if idx(1) == 1
    A = cat(2, 0, A);
end
if idx(end) == 1
    B = cat(2, B, length(idx));
end
Z = find(B-A > thr);
idx_0 = zeros(1, length(idx));
idx_0(A(Z)+1 : B(Z)) = 1;
idx_0 = logical(idx_0);
end
