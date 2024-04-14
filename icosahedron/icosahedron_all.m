% Calculate r values for various N
clear;
close all
tolerance = 1e-8;  % Tolerance for degenerate energy levels
corner_potential = 0;
energy_cut_factor = 3;

% Find max N
upper_lim = 30000; % Max number of sites allowed in laptop memory
max_N = floor(sqrt(upper_lim / 5) - 1);
fprintf(['Max N = ' num2str(max_N) '\n'])
N_nums = primes(max_N);
N_nums = N_nums(5:end);
%N_nums = [23];

% Make directory
folderName = ['icosahedron_plots_ecut=' num2str(energy_cut_factor)];
mkdir(folderName);

r_array = zeros(size(N_nums, 1), 11);
size_array = zeros(size(N_nums, 1), 11);
solvable_prop = zeros(size(N_nums, 1), 2);

index_n = 1;
index_r = 1;
index_size = 1;


for n = 1:size(N_nums, 2)
    N = N_nums(n);
    total_sites = 5 * (N + 1)^2;

    %% Diagonalize H and sigma_d matrices
    tic;
    fprintf(['\n    N = ' num2str(N) '\n'])

    % Generate Hamiltonian matrix
    L = 1/N; % Spacing
    H = zeros(total_sites, total_sites); % Hamiltonian matrix
    
    for i = 1:total_sites
        %fprintf([num2str(i) ' '])
        H(i, i) = 4/L^2;
        %fprintf([num2str(right(i, N)) ' '])
        H(i, right(i, N)) = -4/(3 * L^2);
        %fprintf([num2str(left(i, N)) ' '])
        H(i, left(i, N)) = -4/(3 * L^2);
        if (has_upper(i, N))
            %fprintf('upper ')
        end
        if (has_upper(i, N))
            %fprintf([num2str(upper(i, N)) ' \n'])
            H(i, upper(i, N)) = -4/(3 * L^2);
        else
            %fprintf([num2str(i) ' ' num2str(lower(i, N)) ' \n'])
            H(i, lower(i, N)) = -4/(3 * L^2);
        end
    end

    % Add corner potentials
    faces = [6, 8, 10, 12, 14, 16, 17, 18, 19, 20];
    for i = 1:size(faces, 2)
        face = faces(i);

        lower_left_corner = index_from_coord(1, 1, face, N); 
        lower_right_corner = index_from_coord(N, 1, face, N); 
        top_corner = index_from_coord(1, (N + 1) / 2, face, N);

        %fprintf([num2str(lower_left_corner) ' ' num2str(lower_right_corner) ' ' num2str(top_corner) '\n'])

        H(lower_left_corner, lower_left_corner) = H(lower_left_corner, lower_left_corner) + corner_potential; 
        H(lower_right_corner, lower_right_corner) = H(lower_right_corner, lower_right_corner)+ corner_potential;
        H(top_corner, top_corner) = H(top_corner, top_corner) + corner_potential;
    end

    % Face 2
    faces = [1, 2, 3, 4, 5, 7, 9, 11, 13, 15];
    for i = 1:size(faces, 2)
        face = faces(i);

        bottom_corner = index_from_coord(1, 1, face, N);
        top_right_corner = index_from_coord(N, (N + 1) / 2, face, N);
        top_left_corner = index_from_coord(1, (N + 1) / 2, face, N);

        H(bottom_corner, bottom_corner) = H(bottom_corner, bottom_corner) + corner_potential; 
        H(top_right_corner, top_right_corner) = H(top_right_corner, top_right_corner)  + corner_potential;
        H(top_left_corner, top_left_corner) = H(top_left_corner, top_left_corner) + corner_potential;
    end
    
    % Generate matrix for sigma_d symmetry 
    sigma_d = zeros(total_sites, total_sites);
    for i = 1:total_sites
        sigma_d(i, sigma_d_legit_index(i, N)) = 1;
    end
    
    % Calculate eigenvectors and eigenvalues
    [eigenvectors, eigenvalues] = eig(H + ((1 + sqrt(5))/2)*sigma_d);
    fprintf('evals and evecs done \n')
    toc
    
    %% Generate Symmetry Matrices
    
    % Generate matrix for C3
    C3 = zeros(total_sites, total_sites);
    
    for i = 1: (total_sites)
        C3(i, c3_index(i, N)) = 1;
    end
    
    % Generate matrix for c5
    C5 = zeros(total_sites, total_sites);
    for i = 1: (total_sites)
        C5(i, c5_index(i, N)) = 1;
    end

    % Generate matrix for inv
    inv = zeros(total_sites, total_sites);
    for i = 1: (total_sites)
        inv(i, inv_index(i, N)) = 1;
    end
    
    
    %% Calculate characters of symmetries
    tic;
    % Preallocate data array
    sigma_d_evals = zeros(total_sites, 1);
    c3_evals = zeros(total_sites, 1);
    c5_evals = zeros(total_sites, 1);
    inv_evals = zeros(total_sites, 1);
    
    % Subtract off sigma2 eigenvalues from eigenvalues 
    for i = 1:(total_sites)
        sigma_d_evals(i) = (eigenvectors(:, i).' * sigma_d * eigenvectors(:, i));
        eigenvalues(i, i) = eigenvalues(i, i) - ((1 + sqrt(5))/2) * sigma_d_evals(i);
    end
    fprintf('sigma_d done \n')

    % Determine 1/3rd of max energy
    max_energy = eigenvalues(total_sites, total_sites);
    max_energy_index = 0;

    for i = 1:total_sites
        if (eigenvalues(i, i) <= max_energy / energy_cut_factor)
            max_energy_index = max_energy_index + 1;
        end
    end

    % Compute traces of other symmetries
    for i = 1:max_energy_index
        c5_evals(i) = (eigenvectors(:, i).' * C5 * eigenvectors(:, i));
    end
    fprintf('inv done \n')
    
    for i = 1:max_energy_index
        c3_evals(i) = (eigenvectors(:, i).' * C3 * eigenvectors(:, i));
    end
    fprintf('c3 done \n')

    for i = 1:max_energy_index
        inv_evals(i) = (eigenvectors(:, i).' * inv * eigenvectors(:, i));
    end
    fprintf('inv done \n')
    
    toc;
    
    %% Analyze Degeneracies
    % Reorder energies and evecs
    energy_levels = zeros(max_energy_index, 7);
    [energies, id] = sort(diag(eigenvalues));
    eigenvectors = eigenvectors(:, id);
    sigma_d_evals_sorted = sigma_d_evals(id);
    c3_evals_sorted = c3_evals(id);
    c5_evals_sorted = c5_evals(id);
    inv_evals_sorted = inv_evals(id);
    
    % Analyze for degeneracies
    energy = energies(1);
    degeneracy = 1;
    index = 1;
    trace_sigma2 = sigma_d_evals_sorted(1);
    trace_c3 = c3_evals_sorted(1);
    trace_c5 = c5_evals_sorted(1);
    trace_inv = inv_evals_sorted(1);
    
    for i = 2:max_energy_index
        if (abs(energies(i) - energy) < tolerance) % Next energy is degenerate
            degeneracy = degeneracy + 1;
            trace_sigma2 = trace_sigma2 + sigma_d_evals_sorted(i);
            trace_c3 = trace_c3 + c3_evals_sorted(i);
            trace_c5 = trace_c5 + c5_evals_sorted(i);
            trace_inv = trace_inv + inv_evals_sorted(i);
        else % Next energy is new
            % Record stats for previous energy level
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_sigma2;
            energy_levels(index, 4) = trace_c3;
            energy_levels(index, 5) = trace_c5;
            energy_levels(index, 6) = trace_inv;
            energy_levels(index, 7) = i - 1;
    
            % Reset variables for new energy level
            energy = energies(i);
            degeneracy = 1;
            index = index + 1;
            trace_sigma2 = sigma_d_evals_sorted(i);
            trace_c3 = c3_evals_sorted(i);
            trace_c5 = c5_evals_sorted(i);
            trace_inv = inv_evals_sorted(i);
        end
    
        % Record energy level if we reach the end
        if (i == max_energy_index)
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_sigma2;
            energy_levels(index, 4) = trace_c3;
            energy_levels(index, 5) = trace_c5;
            energy_levels(index, 6) = trace_inv;
            energy_levels(index, 7) = i;
        end
    end
    
    energy_levels = energy_levels(1:index, :);
    
    %% Separate the representation classes
    % Allocate space for each irrep
    elevels_ag = zeros(max_energy_index, 7);
    elevels_t1g = zeros(max_energy_index, 7);
    elevels_t2g = zeros(max_energy_index, 7);
    elevels_gg = zeros(max_energy_index, 7);
    elevels_hg = zeros(max_energy_index, 7);
    elevels_au = zeros(max_energy_index, 7);
    elevels_t1u = zeros(max_energy_index, 7);
    elevels_t2u = zeros(max_energy_index, 7);
    elevels_gu = zeros(max_energy_index, 7);
    elevels_hu = zeros(max_energy_index, 7);
    
    index_ag = 1;
    index_t1g = 1;
    index_t2g = 1;
    index_gg = 1;
    index_hg = 1;
    index_au = 1;
    index_t1u = 1;
    index_t2u = 1;
    index_gu = 1;
    index_hu = 1;
    
    % Round the energy_levels characters
    energy_levels_rounded = zeros(size(energy_levels, 1), 6);
    energy_levels_rounded(:, 2:4) = round(energy_levels(:, 2:4));
    energy_levels_rounded(:, 6) = round(energy_levels(:, 6));
    energy_levels_rounded(:, 1) = energy_levels(:, 1);
    energy_levels_rounded(:, 7) = energy_levels(:, 7);

    % Round the C5 traces
    for i = 1:size(energy_levels, 1)
        if (abs(energy_levels(i, 5) + 2 * cos(2 * pi / 5)) < 0.01)
            energy_levels_rounded(i, 5) = - 2 * cos(2 * pi / 5);
        elseif (abs(energy_levels(i, 5) + 2 * cos(4 * pi / 5)) < 0.01)
            energy_levels_rounded(i, 5) = - 2 * cos(4 * pi / 5);
        else
            energy_levels_rounded(i, 5) = round(energy_levels(i, 5));
        end
    end
    
    % Fill in elevel info in each irrep
    for i = 1:size(energy_levels_rounded, 1)
        trace_e = energy_levels_rounded(i, 2);
        trace_sigma_d = energy_levels_rounded(i, 3);
        trace_c3 = energy_levels_rounded(i, 4);
        trace_c5 = energy_levels_rounded(i, 5);
        trace_inv = energy_levels_rounded(i, 6);
    
        traces = [trace_e, trace_sigma_d, trace_c3, trace_c5, trace_inv];
    
        if (isequal(traces, [1, 1, 1, 1, 1])) % Ag
            elevels_ag(index_ag, :) = energy_levels_rounded(i, :);
            index_ag = index_ag + 1;
        elseif (isequal(traces, [3, -1, 0, -2 * cos(4 * pi / 5), 3])) % T1g
            elevels_t1g(index_t1g, :) = energy_levels_rounded(i, :);
            index_t1g = index_t1g + 1;
        elseif (isequal(traces, [3, -1, 0, -2 * cos(2 * pi / 5), 3])) % T2g
            elevels_t2g(index_t2g, :) = energy_levels_rounded(i, :);
            index_t2g = index_t2g + 1;
        elseif (isequal(traces, [4, 0, 1, -1, 4])) % Gg
            elevels_gg(index_gg, :) = energy_levels_rounded(i, :);
            index_gg = index_gg + 1;
        elseif (isequal(traces, [5, 1, -1, 0, 5])) % Hg
            elevels_hg(index_hg, :) = energy_levels_rounded(i, :);
            index_hg = index_hg + 1;
        elseif (isequal(traces, [1, -1, 1, 1, -1])) % Au
            elevels_au(index_au, :) = energy_levels_rounded(i, :);
            index_au = index_au + 1;
        elseif (isequal(traces, [3, 1, 0, -2 * cos(4 * pi / 5), -3])) % T1u
            elevels_t1u(index_t1u, :) = energy_levels_rounded(i, :);
            index_t1u = index_t1u + 1;
        elseif (isequal(traces, [3, 1, 0, -2 * cos(2 * pi / 5), -3])) % T2u
            elevels_t2u(index_t2u, :) = energy_levels_rounded(i, :);
            index_t2u = index_t2u + 1;
        elseif (isequal(traces, [4, 0, 1, -1, -4])) % Gu
            elevels_gu(index_gu, :) = energy_levels_rounded(i, :);
            index_gu = index_gu + 1;
        elseif (isequal(traces, [5, -1, -1, 0, -5])) % Hu
            elevels_hu(index_hu, :) = energy_levels_rounded(i, :);
            index_hu = index_hu + 1;
        elseif (isequal(traces, [2, 0, 2, 2, 0])) % Accidental degeneracy Ag + Au
            elevels_ag(index_ag, :) = energy_levels_rounded(i, :);
            index_ag = index_ag + 1;
            elevels_au(index_au, :) = energy_levels_rounded(i, :);
            index_au = index_au + 1;
        else 
            fprintf([num2str(energy_levels_rounded(i, 1)) ' ' num2str(energy_levels_rounded(i, 2)) ' ' ...
                num2str(energy_levels_rounded(i, 3)) ' ' num2str(energy_levels_rounded(i, 4)) ' ' ...
                num2str(energy_levels_rounded(i, 5)) ' ' num2str(energy_levels_rounded(i, 6)) '\n'])
        end
    end
    
    % Remove extra rows of zeros
    if (index_ag > 1)
        elevels_ag = elevels_ag(1:(index_ag - 1), :);
    end 
    
    if (index_t1g > 1)
        elevels_t1g = elevels_t1g(1:(index_t1g - 1), :);
    end
    
    if (index_t2g > 1)
        elevels_t2g = elevels_t2g(1:(index_t2g - 1), :);
    end
    
    if (index_gg > 1)
        elevels_gg = elevels_gg(1:(index_gg - 1), :);
    end
    
    if (index_hg > 1)
        elevels_hg = elevels_hg(1:(index_hg - 1), :);
    end

    if (index_au > 1)
        elevels_au = elevels_au(1:(index_au - 1), :);
    end 
    
    if (index_t1u > 1)
        elevels_t1u = elevels_t1u(1:(index_t1u - 1), :);
    end
    
    if (index_t2u > 1)
        elevels_t2u = elevels_t2u(1:(index_t2u - 1), :);
    end
    
    if (index_gu > 1)
        elevels_gu = elevels_gu(1:(index_gu - 1), :);
    end
    
    if (index_hu > 1)
        elevels_hu = elevels_hu(1:(index_hu - 1), :);
    end
    
    %% Find energy level spacings
    % Make vector of energy level spacings
    spacings_ag = zeros(size(elevels_ag, 1) - 1, 1);
    spacings_t1g = zeros(size(elevels_t1g, 1) - 1, 1);
    spacings_t2g = zeros(size(elevels_t2g, 1) - 1, 1);
    spacings_gg = zeros(size(elevels_gg, 1) - 1, 1);
    spacings_hg = zeros(size(elevels_hg, 1) - 1, 1);
    spacings_au = zeros(size(elevels_au, 1) - 1, 1);
    spacings_t1u = zeros(size(elevels_t1u, 1) - 1, 1);
    spacings_t2u = zeros(size(elevels_t2u, 1) - 1, 1);
    spacings_gu = zeros(size(elevels_gu, 1) - 1, 1);
    spacings_hu = zeros(size(elevels_hu, 1) - 1, 1);
    
    for i = 1:(size(elevels_ag, 1) - 1)
        spacings_ag(i) = abs(elevels_ag(i, 1) - elevels_ag(i+1, 1));
    end
    
    for i = 1:(size(elevels_t1g, 1) - 1)
        spacings_t1g(i) = abs(elevels_t1g(i, 1) - elevels_t1g(i+1, 1));
    end
    
    for i = 1:(size(elevels_t2g, 1) - 1)
        spacings_t2g(i) = abs(elevels_t2g(i, 1) - elevels_t2g(i+1, 1));
    end
    
    for i = 1:(size(elevels_gg, 1) - 1)
        spacings_gg(i) = abs(elevels_gg(i, 1) - elevels_gg(i+1, 1));
    end
    
    for i = 1:(size(elevels_hg, 1) - 1)
        spacings_hg(i) = abs(elevels_hg(i, 1) - elevels_hg(i+1, 1));
    end
    
    for i = 1:(size(elevels_au, 1) - 1)
        spacings_au(i) = abs(elevels_au(i, 1) - elevels_au(i+1, 1));
    end
    
    for i = 1:(size(elevels_t1u, 1) - 1)
        spacings_t1u(i) = abs(elevels_t1u(i, 1) - elevels_t1u(i+1, 1));
    end
    
    for i = 1:(size(elevels_t2u, 1) - 1)
        spacings_t2u(i) = abs(elevels_t2u(i, 1) - elevels_t2u(i+1, 1));
    end
    
    for i = 1:(size(elevels_gu, 1) - 1)
        spacings_gu(i) = abs(elevels_gu(i, 1) - elevels_gu(i+1, 1));
    end
    
    for i = 1:(size(elevels_hu, 1) - 1)
        spacings_hu(i) = abs(elevels_hu(i, 1) - elevels_hu(i+1, 1));
    end

    % Compute level spacing ratios
    r_ag = LSR(spacings_ag);
    r_t1g = LSR(spacings_t1g);
    r_t2g = LSR(spacings_t2g);
    r_gg = LSR(spacings_gg);
    r_hg = LSR(spacings_hg);
    r_au = LSR(spacings_au);
    r_t1u = LSR(spacings_t1u);
    r_t2u = LSR(spacings_t2u);
    r_gu = LSR(spacings_gu);
    r_hu = LSR(spacings_hu);

    % Fill in r array values
    r_array(index_r, :) = [N, r_ag, r_t1g, r_t2g, r_gg, r_hg, r_au, r_t1u, r_t2u, r_gu, r_hu];
   
    % Fill in size of irreps
    size_array(index_size, :) = [N, size(elevels_ag, 1), size(elevels_t1g, 1), size(elevels_t2g, 1), ...
        size(elevels_gg, 1), size(elevels_hg, 1), size(elevels_au, 1), size(elevels_t1u, 1), size(elevels_t2u, 1), ...
        size(elevels_gu, 1), size(elevels_hu, 1)];

    % Find proportion of solvable states
    solvable_prop(index_size, :) = [size(elevels_ag, 1), size(elevels_t1g, 1)] / (total_sites);

    %% Plot the level spacings histogram and show r values
    % Plot histograms
    figure(index_n)
    tiledlayout(2, 5,'TileSpacing', 'tight','Padding','Tight')
    bin_factor = 5;
    
    nexttile
    histogram(spacings_ag, ceil(size(elevels_ag, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Ag')
    subtitle(['r = ' num2str(r_ag) '; total = ' num2str(size_array(index_size, 2))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t1g, ceil(size(elevels_t1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('T1g')
    subtitle(['r = ' num2str(r_t1g) '; total = ' num2str(size_array(index_size, 3))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t2g, ceil(size(elevels_t2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('T2g')
    subtitle(['r = ' num2str(r_t2g) '; total = ' num2str(size_array(index_size, 4))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_gg, ceil(size(elevels_gg, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Gg')
    subtitle(['r = ' num2str(r_gg) '; total = ' num2str(size_array(index_size, 5))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_hg, ceil(size(elevels_hg, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Hg')
    subtitle(['r = ' num2str(r_hg) '; total = ' num2str(size_array(index_size, 6))])
    xlabel('Energy spacing')
    ylabel('Counts')

    nexttile
    histogram(spacings_au, ceil(size(elevels_au, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Au')
    subtitle(['r = ' num2str(r_au) '; total = ' num2str(size_array(index_size, 7))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t1u, ceil(size(elevels_t1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('T1u')
    subtitle(['r = ' num2str(r_t1u) '; total = ' num2str(size_array(index_size, 8))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t2u, ceil(size(elevels_t2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('T2u')
    subtitle(['r = ' num2str(r_t2u) '; total = ' num2str(size_array(index_size, 9))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_gu, ceil(size(elevels_gu, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Gu')
    subtitle(['r = ' num2str(r_gu) '; total = ' num2str(size_array(index_size, 10))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_hg, ceil(size(elevels_hg, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Hg')
    subtitle(['r = ' num2str(r_hg) '; total = ' num2str(size_array(index_size, 11))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    
    set(figure(index_n),'position',[0,100,3000,400])
    saveas(gcf, [folderName '/N=' num2str(N) '_elevel_hist.jpeg']);
    
    % Increment indices
    index_n = index_n + 1;
    index_size = index_size + 1;
    index_r = index_r + 1;
end

%% Plot r as a function of size in each irrep
figure(index_n)
tiledlayout(2, 5, 'TileSpacing', 'tight', 'Padding', 'Tight')
ylow = 0.3;
yhigh = 0.7;

nexttile
plot(N_nums.', r_array(:, 2).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Ag')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 3).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('T1g')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 4).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('T2g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 5).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Gg')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 6).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Hg')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

nexttile
plot(N_nums.', r_array(:, 7).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Au')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 8).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('T1u')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 9).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('T2u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 10).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Gu')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 11).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Hu')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow, yhigh])

set(figure(index_n),'position',[0,100,3000,400])
saveas(gcf, [folderName '/LSR_plot.jpeg']);

%% Plot proportions of solvable states
figure(index_n + 1)
plot(N_nums, sum(solvable_prop, 2), 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Proportion of states which are solvable')
xlabel('N')
yline(1/12, '--', '1/12')
ylim([0.0, 0.1])
ylabel('Proportion')
saveas(gcf, [folderName '/solv_prop.jpeg']);

%% Save variables
save([folderName '/r_array.mat'],'r_array','-v7.3')
save([folderName '/size_array.mat'],'size_array','-v7.3')
save([folderName '/elevels_ag.mat'],'elevels_ag','-v7.3')
save([folderName '/elevels_t1g.mat'],'elevels_t1g','-v7.3')
save([folderName '/elevels_t2g.mat'],'elevels_t2g','-v7.3')
save([folderName '/elevels_gg.mat'],'elevels_gg','-v7.3')
save([folderName '/elevels_hg.mat'],'elevels_hg','-v7.3')
save([folderName '/elevels_au.mat'],'elevels_au','-v7.3')
save([folderName '/elevels_t1u.mat'],'elevels_t1u','-v7.3')
save([folderName '/elevels_t2u.mat'],'elevels_t2u','-v7.3')
save([folderName '/elevels_gu.mat'],'elevels_gu','-v7.3')
save([folderName '/elevels_hu.mat'],'elevels_hu','-v7.3')


%% Functions

% Gives width of row within face CHECKED
function w = row_width(index, N) 
    face = face_from_index(index, N);
    y = yfacecoord(index, N);
    if ((face == 1) || (face == 2) || (face == 3) || (face == 4) || (face == 5) || (face == 7)...
            || (face == 9) || (face == 11) || (face == 13) || (face == 15))
        w = 1 + 2 * (y - 1);
    else
        w = N - 2 * (y - 1);
    end
end

% Gives y coordinate of index within each face CHECKED
function y = yfacecoord(index, N) 
    sites_per_face = ((N + 1) / 2)^2;
    face = face_from_index(index, N);
    if (index <= 5 * sites_per_face) % bottom row
        index = mod(index - 1, sites_per_face) + 1;
        y = ceil(sqrt(index));
    elseif (index > 15 * sites_per_face)
        total_sites = 5 * (N + 1)^2 - (20 - face) * sites_per_face;
        y = (N + 1) / 2 + 1 - ceil(sqrt(total_sites + 1 - index));
    else
        % remove bottom row of triangle
        index = index - 5 * sites_per_face;
        row_length = 5 * (N + 1);

        y = floor((index - 1) / row_length) + 1;
    end
end

% Gives x coordinate of index within each face CHECKED
function x = xfacecoord(index, N) 
    face = face_from_index(index, N);
    y_coord = yfacecoord(index, N);
    sites_per_face = ((N + 1) / 2)^2;
    row_w = row_width(index, N);

    if (face <= 5) % bottom row
        index = mod(index - 1, sites_per_face) + 1;
        x = index - (y_coord - 1)^2; 
    elseif (face > 15) %top row
        total_sites = 5 * (N + 1)^2 - (20 - face) * sites_per_face;
        x = row_w - (total_sites + 1 - index - ((N + 1) / 2 - y_coord)^2) + 1;
    elseif (face == 6)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1);
    elseif (face == 8)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - (N + 1);
    elseif (face == 10)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - 2 * (N + 1);
    elseif (face == 12)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - 3 * (N + 1);
    elseif (face == 14)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - 4 * (N + 1);
    elseif (face == 7)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - (N + 1 - row_w);
    elseif (face == 9)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - (N + 1 - row_w) - (N + 1);
    elseif (face == 11)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - (N + 1 - row_w) - 2 * (N + 1);
    elseif (face == 13)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - (N + 1 - row_w) - 3 * (N + 1);
    elseif (face == 15)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - (N + 1 - row_w) - 4 * (N + 1);
    end
end

% Gives face containing given index CHECKED
function f = face_from_index(index, N) 
    sites_per_face = ((N + 1) / 2)^2;
    if (index <= 5 * sites_per_face) % bottom row of triangles
        f = ceil(index / sites_per_face);
    elseif (index > 15 * sites_per_face) % top row of triangles
        f = ceil(index / sites_per_face);
    else % in middle row of triangles
        % remove bottom row of triangle
        index = index - 5 * sites_per_face;
        row_length = 5 * (N + 1);

        face_pair = mod(ceil(index / (N + 1)) - 1, 5) + 1;
        row_num = ceil(index / row_length);

        short_width = 1 + 2 * (row_num - 1);
        long_width = (N + 1) - short_width;

        index = mod(index - 1, short_width + long_width) + 1;
        if (index <= long_width)
            f = 2 * face_pair + 4; 
        else
            f = 2 * face_pair + 5;
        end
    end
end

% From face coordinate and face number give index CHECKED
function i = index_from_coord(x, y, face, N)
    sites_per_face = ((N + 1) / 2)^2;

    if (face <= 5)
        i = sites_per_face * (face - 1) + (y - 1)^2 + x;
    elseif (face == 6)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + x;
    elseif (face == 7)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + (N + 1 - (1 + 2 * (y - 1))) + x;
    elseif (face == 8)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + (N + 1) + x;
    elseif (face == 9)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + (N + 1) + (N + 1 - (1 + 2 * (y - 1))) + x;
    elseif (face == 10)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + 2 * (N + 1) + x;
    elseif (face == 11)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + 2 * (N + 1) + (N + 1 - (1 + 2 * (y - 1))) + x;
    elseif (face == 12)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + 3 * (N + 1) + x;
    elseif (face == 13)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + 3 * (N + 1) + (N + 1 - (1 + 2 * (y - 1))) + x;
    elseif (face == 14)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + 4 * (N + 1) + x;
    elseif (face == 15)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + 4 * (N + 1) + (N + 1 - (1 + 2 * (y - 1))) + x;
    elseif (face > 15)
        total_sites = 5 * (N + 1)^2 - (20 - face) * sites_per_face;
        i = total_sites - ((N + 1) / 2 + 1 - y)^2 + x;
    end
end

% Given index on discretized octahedron returns new index after C3 symmetry
% CHECKED
function i = c3_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);
    if (face == 1) 
        new_face = 13;
        new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;
        new_y = (N + 1) / 2 - floor((x_coord - 1) / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 12;
        new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 11;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 - y_coord + 1 + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 18;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = 1 + floor((x_coord - 1) / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 5) 
        new_face = 19;
        new_x = x_coord;
        new_y = y_coord - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 6) 
        new_face = 14;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 - y_coord + 1 - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 7) 
        new_face = 5;
        new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;
        new_y = (N + 1) / 2 - floor((x_coord - 1) / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 8) 
        new_face = 4;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 9) 
        new_face = 3;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 - y_coord + 1 + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 10) 
        new_face = 10;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = 1 + floor((x_coord - 1) / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 11) 
        new_face = 9;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 - y_coord + 1 + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 12) 
        new_face = 17;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = 1 + floor((x_coord - 1) / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 13) 
        new_face = 16;
        new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 14) 
        new_face = 20;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 - y_coord + 1 - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 15) 
        new_face = 15;
        new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;
        new_y = (N + 1) / 2 - floor((x_coord - 1) / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 16) 
        new_face = 1;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 - floor((x_coord - 1) / 2); 
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 17) 
        new_face = 2;
        new_x = x_coord;
        new_y = y_coord + floor(x_coord / 2);
        %fprintf(['x=' num2str(new_x) ' y=' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 18) 
        new_face = 8;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = 1 + floor((x_coord - 1) / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 19) 
        new_face = 7;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 20) 
        new_face = 6;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 - y_coord + 1 - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized octahedron returns new index after sigma_d symmetry
% CHECKED
function i = sigma_d_legit_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);
    if (face == 1) 
        new_face = 1;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 5;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 4;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 3;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 5) 
        new_face = 2;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 6) 
        new_face = 6;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 7) 
        new_face = 15;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 8) 
        new_face = 14;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 9) 
        new_face = 13;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 10) 
        new_face = 12;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 11) 
        new_face = 11;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 12) 
        new_face = 10;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 13) 
        new_face = 9;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 14) 
        new_face = 8;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 15) 
        new_face = 7;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 16) 
        new_face = 20;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 17) 
        new_face = 19;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 18) 
        new_face = 18;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 19) 
        new_face = 17;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 20) 
        new_face = 16;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized octahedron returns new index after C5 symmetry
% CHECKED
function i = c5_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    if (face == 1) 
        new_face = 14;
        new_x = x_coord;
        new_y = y_coord - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 5;
        new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;
        new_y = (N + 1) / 2 - floor((x_coord - 1) / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 4;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 - y_coord + 1 + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 12;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = 1 + floor((x_coord - 1) / 2);
        %fprintf(['x=' num2str(new_x) ' y=' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 5) 
        new_face = 13;
        new_x = x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 6) 
        new_face = 15;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 - floor((x_coord - 1) / 2); 
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 7) 
        new_face = 6;
        new_x = x_coord;
        new_y = y_coord - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 8) 
        new_face = 1;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 - floor((x_coord - 1) / 2); 
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 9) 
        new_face = 2;
        new_x = x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 10) 
        new_face = 3;
        new_x = x_coord;
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 11) 
        new_face = 10;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = 1 + floor((x_coord - 1) / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 12) 
        new_face = 11;
        new_x = x_coord;
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 13) 
        new_face = 18;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = 1 + floor((x_coord - 1) / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 14) 
        new_face = 19;
        new_x = x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 15) 
        new_face = 20;
        new_x = x_coord;
        new_y = y_coord - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 16) 
        new_face = 7;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 - floor((x_coord - 1) / 2); 
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 17) 
        new_face = 8;
        new_x = x_coord;
        new_y = y_coord;
        %fprintf(['x=' num2str(new_x) ' y=' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 18) 
        new_face = 9;
        new_x = x_coord;
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 19) 
        new_face = 17;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = 1 + floor((x_coord - 1) / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 20) 
        new_face = 16;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 - y_coord + 1 - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized octahedron returns new index after inv
% symmetry CHECKED
function i = inv_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);
    if (face == 1) 
        new_face = 18;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 19;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 20;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 16;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 5) 
        new_face = 17;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 6) 
        new_face = 11;
        new_x = x_coord;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 7) 
        new_face = 12;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 8) 
        new_face = 13;
        new_x = x_coord;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 9) 
        new_face = 14;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 10) 
        new_face = 15;
        new_x = x_coord;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 11) 
        new_face = 6;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 12) 
        new_face = 7;
        new_x = x_coord;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 13) 
        new_face = 8;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 14) 
        new_face = 9;
        new_x = x_coord;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 15) 
        new_face = 10;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 16) 
        new_face = 4;
        new_x = x_coord;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 17) 
        new_face = 5;
        new_x = x_coord;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 18) 
        new_face = 1;
        new_x = x_coord;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 19) 
        new_face = 2;
        new_x = x_coord;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 20) 
        new_face = 3;
        new_x = x_coord;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Functions to find indices of H matrix 
function y = upper(index, N) % CHECKED
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face <= 5 && y_coord < (N + 1) / 2) % in faces 1-5 but not top row
        y = index_from_coord(x_coord + 1, y_coord + 1, face, N);
    elseif (face == 1 && y_coord == (N + 1) / 2) % in face 1 top row
        y = index_from_coord(x_coord, 1, 6, N);
    elseif (face == 2 && y_coord == (N + 1) / 2) % in face 2 top row
        y = index_from_coord(x_coord, 1, 8, N);
    elseif (face == 3 && y_coord == (N + 1) / 2) % in face 3 top row
        y = index_from_coord(x_coord, 1, 10, N);
    elseif (face == 4 && y_coord == (N + 1) / 2) % in face 4 top row
        y = index_from_coord(x_coord, 1, 12, N);
    elseif (face == 5 && y_coord == (N + 1) / 2) % in face 5 top row
        y = index_from_coord(x_coord, 1, 14, N);
    elseif (face > 15) % faces 15-20
        y = index_from_coord(x_coord - 1, y_coord + 1, face, N);
    elseif (face == 7 && y_coord == (N + 1) / 2) % face 7 top row
        y = index_from_coord(x_coord, 1, 16, N);
    elseif (face == 9 && y_coord == (N + 1) / 2) % face 9 top row
        y = index_from_coord(x_coord, 1, 17, N);
    elseif (face == 11 && y_coord == (N + 1) / 2) % face 11 top row
        y = index_from_coord(x_coord, 1, 18, N);
    elseif (face == 13 && y_coord == (N + 1) / 2) % face 13 top row
        y = index_from_coord(x_coord, 1, 19, N);
    elseif (face == 15 && y_coord == (N + 1) / 2) % face 15 top row
        y = index_from_coord(x_coord, 1, 20, N);
    elseif (mod(face, 2) == 0) % face = 6, 8, 10, 12, 14 
        y = index_from_coord(x_coord - 1, y_coord + 1, face, N);
    else % face = 7, 9, 11, 13, 15 not rop row
        y = index_from_coord(x_coord + 1, y_coord + 1, face, N);
    end
end

function y = lower(index, N) % CHECKED
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face <= 5)
        y = index_from_coord(x_coord - 1, y_coord - 1, face, N);
    elseif (face == 6 && y_coord == 1) % face 6 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 1, N);
    elseif (face == 8 && y_coord == 1) % face 8 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 2, N);
    elseif (face == 10 && y_coord == 1) % face 10 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 3, N);
    elseif (face == 12 && y_coord == 1) % face 12 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 4, N);
    elseif (face == 14 && y_coord == 1) % face 14 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 5, N);
    elseif (face == 16 && y_coord == 1) % face 16 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 7, N);
    elseif (face == 17 && y_coord == 1) % face 17 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 9, N);
    elseif (face == 18 && y_coord == 1) % face 18 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 11, N);
    elseif (face == 19 && y_coord == 1) % face 19 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 13, N);
    elseif (face == 20 && y_coord == 1) % face 20 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 15, N);
    elseif (face > 15)
        y = index_from_coord(x_coord + 1, y_coord - 1, face, N);
    else
        y = index - 5 * (N + 1) + 1;
    end
end

function y = right(index, N) % CHECKED
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);

    if (face == 1 && x_coord == row_w) % face 1 right side
        y = index_from_coord(1, y_coord, 2, N);
    elseif (face == 2 && x_coord == row_w) % face 2 right side
        y = index_from_coord(1, y_coord, 3, N);
    elseif (face == 3 && x_coord == row_w) % face 3 right side
        y = index_from_coord(1, y_coord, 4, N);
    elseif (face == 4 && x_coord == row_w) % face 4 right side
        y = index_from_coord(1, y_coord, 5, N);
    elseif (face == 5 && x_coord == row_w) % face 5 right side
        y = index_from_coord(1, y_coord, 1, N);
    elseif (face == 16 && x_coord == row_w) % face 16 right side
        y = index_from_coord(1, y_coord, 17, N);
    elseif (face == 17 && x_coord == row_w) % face 17 right side
        y = index_from_coord(1, y_coord, 18, N);
    elseif (face == 18 && x_coord == row_w) % face 18 right side
        y = index_from_coord(1, y_coord, 19, N);
    elseif (face == 19 && x_coord == row_w) % face 19 right side
        y = index_from_coord(1, y_coord, 20, N);
    elseif (face == 20 && x_coord == row_w) % face 20 right side
        y = index_from_coord(1, y_coord, 16, N);
    elseif (face == 15 && x_coord == row_w) % face 7 right side
        y = index_from_coord(1, y_coord, 6, N);
    else
        y = index + 1; 
    end
end

function y = left(index, N) % CHECKED
    face = face_from_index(index, N);
    row_w = row_width(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face == 1 && x_coord == 1) % face 1 left side
        y = index_from_coord(row_w, y_coord, 5, N);
    elseif (face == 2 && x_coord == 1) % face 2 left side
        y = index_from_coord(row_w, y_coord, 1, N);
    elseif (face == 3 && x_coord == 1) % face 3 left side
        y = index_from_coord(row_w, y_coord, 2, N);
    elseif (face == 4 && x_coord == 1) % face 4 left side
        y = index_from_coord(row_w, y_coord, 3, N);
    elseif (face == 5 && x_coord == 1) % face 5 left side
        y = index_from_coord(row_w, y_coord, 4, N);
    elseif (face == 16 && x_coord == 1) % face 16 left side
        y = index_from_coord(row_w, y_coord, 20, N);
    elseif (face == 17 && x_coord == 1) % face 17 left side
        y = index_from_coord(row_w, y_coord, 16, N);
    elseif (face == 18 && x_coord == 1) % face 18 left side
        y = index_from_coord(row_w, y_coord, 17, N);
    elseif (face == 19 && x_coord == 1) % face 19 left side
        y = index_from_coord(row_w, y_coord, 18, N);
    elseif (face == 20 && x_coord == 1) % face 20 left side
        y = index_from_coord(row_w, y_coord, 19, N);
    elseif (face == 6 && x_coord == 1)
        y = index_from_coord(N + 1 - row_w, y_coord, 15, N);
    else
        y = index - 1;
    end  
end

% Indicate whether nearest vertical neighbor is above or below CHECKED
function answer = has_upper(index, N)
    y = yfacecoord(index, N);
    face = face_from_index(index, N);
    sites_per_face = ((N + 1) / 2)^2;

    if (face <= 5) % in lower row 
        index = mod(index - 1, sites_per_face) + 1; 
        if (mod(y, 2) - mod(index, 2) == 0)
            answer = true;
        else 
            answer = false;
        end
    elseif (face > 15) % in upper row
        index = index - (face - 1) * sites_per_face;
        if (mod(y, 2) - mod(index, 2) == 0)
            answer = false;
        else 
            answer = true;
        end
    else % in middle row
        index = index - 5 * sites_per_face;
        if (mod(index, 2) == 0)
            answer = true;
        else
            answer = false;
        end
    end
end

% Function to calculate mean level spacing ratio
function r = LSR(spacings)
    r_n = zeros(size(spacings, 1) - 1, 1);
    
    for i = 1:(size(spacings, 1) - 1)
        r_n(i) = min(spacings(i), spacings(i + 1)) ...
            / max(spacings(i), spacings(i + 1));
    end

    r = mean(r_n);
end