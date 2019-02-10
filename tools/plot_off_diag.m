function plot_off_diag(C)

N = size(C,1);

for n=1:N
    C(n,n) = nan;
end

imagesc(C)
