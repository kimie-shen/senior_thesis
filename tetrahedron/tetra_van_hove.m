% Van Hove Singularity

N_nums = [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163];
degeneracies_1 = [9, 15, 18, 24, 27, 33, 42, 45, 54, 60, 63, 69, 78, 87, 90, 99, 105, 108, 117, 123, 132, 144, 150, 153, 159, 162, 168, 189, 195, 204, 207, 222, 225, 234, 243];
degeneracies_2 = [12, 18, 21, 27, 30, 36, 45, 48, 57, 63, 66, 72, 81, 90, 93, 102, 108, 111, 120, 126, 135, 147, 153, 156, 162, 165, 171, 192, 198, 207, 210, 225, 228, 237, 246];
energies_1 = [322.6667, 962.6667, 1410.6667, 2562.6667, 4930.6667, ...
    5890.6667, 9282.6667, 11970.6667, 13442.6667];
energies_2 = [645.3333, 1925.3333, 2821.3333, 5125.3333, 9861.3333, ...
    11781.3333, 18565.3333, 23941.3333, 26885.3333];

figure(1)
plot(N_nums, degeneracies_1, 'Linestyle','-','Marker','.')
hold on
plot(N_nums, degeneracies_2, 'Linestyle','-','Marker','.')
xlabel('$N$', 'Interpreter', 'latex')
ylabel('Degeneracy of Van Hove Singularity')

%figure(2)
%tiledlayout(1, 2,'TileSpacing', 'tight','Padding','Tight')
%nexttile
%plot(N_nums, degeneracies, 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
%xlabel('$N$', 'Interpreter', 'latex')
%ylabel('Degeneracy of Van Hove Singularity')

%nexttile
%plot(energies_1, N_nums,'Linestyle','-','Marker','.')
%hold on 
%plot(energies_2, N_nums,'Linestyle','-','Marker','.')
%ylabel('$N$', 'Interpreter', 'latex')
%xlabel('Energy of Van Hove Singularity')
%legend('1st VHS', '2nd VHS', 'Location', 'best')
