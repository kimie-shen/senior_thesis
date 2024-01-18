figure(index_n + 1)
proportions = sum(solvable_prop, 2);
plot(N_nums(2:end), proportions(2:end), '-.')
title('Proportion of states which are solvable')
xlabel('N')
yline(1/4, '--', '1/4')
ylim([0.248, 0.253])
ylabel('Proportion')