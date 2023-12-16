N_nums = [5, 7, 11, 13, 17, 19, 23, 29, 31, 37];
plot(N_nums, solvable_prop, '-o')
title('Proportion of states which are solvable')
xlabel('N')
yline(1/4, '--', '1/4')
ylim([0.245, 0.255])
ylabel('Proportion')