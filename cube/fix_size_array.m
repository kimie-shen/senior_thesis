n_runs = (N_end - N_start) / N_spacing + 1;
solvable_sizes = zeros(n_runs, 4);

for i = 1:n_runs
    N = 5 * i;
    total_size = N / size_array(i, 1);
    fprintf([num2str(total_size) '\n']);
    fprintf([num2str(6 * N^2) '\n \n ']);
    solvable_sizes(i, 1:2) = size_array(i, 2:3) * total_size / (6 * N^2);
    solvable_sizes(i, 3:4) = size_array(i, 7:8) * total_size / (6 * N^2);
end

figure(1)
bar(solvable_sizes, 'stacked')
hold on
yline(1/12)
title('Proportion of Solvable States')
xlabel('N')
ylabel('Proportion')

figure(2)
plot(sum(solvable_sizes, 2), '-o')
hold on
yline(1/12, '--', '1/12')
ylim([0, 0.15])
title('Proportion of Solvable states')
xlabel('N/5')
ylabel('Proportion')