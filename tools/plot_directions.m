function plot_directions(U, n)

if nargin<2
    n = size(U,2);
end


for i=1:n
    hold on
    plot([0; U(1,i)], [0, U(2,i)], 'r', 'linewidth', 2)
end