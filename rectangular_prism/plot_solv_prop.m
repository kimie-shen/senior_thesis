%% Plot proportions of solvable states

% one-dimensional irreps
prop_1d_irreps = zeros(size(size_array, 1), 1);

for i = 1:size(size_array, 1)
    total = size_array(i, 2) + size_array(i, 3) + size_array(i, 4) + size_array(i, 5) ...
        + size_array(i, 6) + size_array(i, 7) + size_array(i, 8) + size_array(i, 9) ...
        + size_array(i, 10) + size_array(i, 11) + size_array(i, 12) + size_array(i, 13);
    prop_1d_irreps(i) = (size_array(i, 10) + size_array(i, 11) + size_array(i, 12) + size_array(i, 13)) / total;
end

N_nums = size_array(:, 1);

figure(1)
plot(N_nums.', prop_1d_irreps, 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Proportion of solvable states')
xlabel('N')
prop_solvable = 1 / (4 * (l1 * l2 + l1 * l3 + l2 * l3));
yline(prop_solvable, '--', ['1/' num2str(4 * (l1 * l2 + l1 * l3 + l2 * l3))])
prop_range = 0.002;
ylim([prop_solvable - prop_range, prop_solvable + prop_range/2])
ylabel('Proportion')
set(figure(1),'position',[0,100,600,350])
%saveas(gcf, [folderName '/solvable_proportion.jpeg']);