% Van Hove Singularity

N_nums = [11, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83, 103];
degeneracies = [20, 36, 44, 60, 84, 92, 116, 132, 140, 156, 164, 204];
energies_1 = [322.6667, 962.6667, 1410.6667, 2562.6667, 4930.6667, ...
    5890.6667, 9282.6667, 11970.6667, 13442.6667, 16642.6667, 18370.6667, 28290.6667];
energies_2 = [645.3333, 1925.3333, 2821.3333, 5125.3333, 9861.3333, ...
    11781.3333, 18565.3333, 23941.3333, 26885.3333, 33285.3333, 36741.3333, 56581.3333];

figure(1)
plot(N_nums, degeneracies, 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
xlabel('$N$', 'Interpreter', 'latex')
ylabel('Degeneracy of Van Hove Singularity')

figure(2)
tiledlayout(1, 2,'TileSpacing', 'tight','Padding','Tight')
nexttile
plot(N_nums, degeneracies, 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
xlabel('$N$', 'Interpreter', 'latex')
ylabel('Degeneracy of Van Hove Singularity')

nexttile
plot(energies_1, N_nums,'Linestyle','-','Marker','.')
hold on 
plot(energies_2, N_nums,'Linestyle','-','Marker','.')
ylabel('$N$', 'Interpreter', 'latex')
xlabel('Energy of Van Hove Singularity')
legend('1st VHS', '2nd VHS', 'Location', 'best')
