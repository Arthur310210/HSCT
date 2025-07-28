function y = Check_zero(x, Thr)
L = length(x);
y = Inf(L,1);
nonzero_inx =find(abs(x)>Thr);
y(nonzero_inx) = x(nonzero_inx);
end