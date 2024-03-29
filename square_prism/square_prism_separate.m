%% Complete calculation for square prism (no corner potential)
clear;
close all
tolerance = 1e-8;  % Tolerance for degenerate energy levels

% Side length ratios
l1 = 2;
l2 = 1;

% Is l1 even or not?
l1_even = false;
if (mod(l1, 2) == 0)
    l1_even = true;
end

% Is l2 even or not?
l2_even = false;
if (mod(l2, 2) == 0)
    l2_even = true;
end

% Find max N
upper_lim = 30000; % Max number of sites allowed in laptop memory
max_N = floor(sqrt(upper_lim / (2 * (l1 * l1) + 4 * (l1 * l2))));
fprintf(['Max N = ' num2str(max_N) '\n'])
N_nums = primes(max_N);
N_nums = N_nums(4:end);

% Make directory
folderName = ['square_prism_l1=' num2str(l1) '_l2=' num2str(l2) '_separated'];
mkdir(folderName);

% Set indices
index_n = 1;
index_r = 1;
index_size = 1;

% Allocate variable space
r_array = zeros(size(N_nums, 1), 15);
size_array = zeros(size(N_nums, 1), 15);
solvable_prop = zeros(size(N_nums, 1), 4);

%% Start loop
for n = 1:size(N_nums, 2)

    N = N_nums(n);
    fprintf(['\n    N = ' num2str(N) '\n'])
    total_sites = (2 * l1^2 + 4 * l1 * l2) * N^2; % Total number of sites in unfolded polyhedron

    %% Diagonalize H and sigma_d matrices
    tic;

    % Generate Hamiltonian matrix
    L = 1/N; % Spacing
    H = zeros(total_sites, total_sites); % Hamiltonian matrix
    
    for i = 1:total_sites
        H(i, i) = 4/L^2;
        H(i, upper(i, N, l1, l2, total_sites)) = -1/L^2;
        H(i, lower(i, N, l1, l2, total_sites)) = -1/L^2;
        H(i, right(i, N, l1, l2, total_sites)) = -1/L^2;
        H(i, left(i, N, l1, l2, total_sites)) = -1/L^2;
    end
    
    % Generate matrix for sigma_d symmetry 
    sigma_d = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        %fprintf([num2str(i) ' ' num2str(sigma_d_index(i, N, l1, l2, total_sites)) '\n'])
        sigma_d(i, sigma_d_index(i, N, l1, l2, total_sites)) = 1;
    end
    
    %% Calculate eigenvectors and eigenvalues
    [eigenvectors, eigenvalues] = eig(H + ((1 + sqrt(5))/2) * sigma_d);
    fprintf('evals and evecs done \n')
    toc
    
    %% Generate Symmetry Matrices
    
    % Generate matrix for C2_oo
    C2_oo = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        C2_oo(i, c2_oo_index(i, N, l1, l2, total_sites)) = 1;
    end

    % Generate matrix for i
    inv = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        inv(i, inv_index(i, N, l1, l2, total_sites)) = 1;
    end

    % Generate matrix for C2_o
    C2_o = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        C2_o(i, c2_o_index(i, N, l1, l2, total_sites)) = 1;
    end
    
    %% Calculate characters of symmetries
    tic;
    % Preallocate data array
    sigma_d_evals = zeros(total_sites, 1);
    c2_oo_evals = zeros(total_sites, 1);
    inv_evals = zeros(total_sites, 1);
    c2_o_evals = zeros(total_sites, 1);
    
    % Find inv eigenvalues
    for i = 1:(total_sites)
        inv_evals(i) = (eigenvectors(:, i).' * inv * eigenvectors(:, i));
    end
    fprintf('inv done \n')
    
    % Find energy and sigma_d eigenvluaes
    for i = 1:(total_sites)
        sigma_d_evals(i) = (eigenvectors(:, i).' * sigma_d * eigenvectors(:, i));
        eigenvalues(i, i) = eigenvalues(i, i) - ((1 + sqrt(5))/2) * sigma_d_evals(i);
    end
    fprintf('sigma_d done \n')
    
    % Find c2_oo eigenvalues
    for i = 1:(total_sites)
        c2_oo_evals(i) = (eigenvectors(:, i).' * C2_oo * eigenvectors(:, i));
    end
    fprintf('c2_oo done \n')
    
    % Find c2_o eigenvalues
    for i = 1:(total_sites)
        c2_o_evals(i) = (eigenvectors(:, i).' * C2_o * eigenvectors(:, i));
    end
    fprintf('c2_o done \n')
    toc;
    
    %% Analyze Degeneracies
    % Reorder energies and evecs
    energy_levels = zeros(total_sites, 7);
    [energies, id] = sort(diag(eigenvalues));
    eigenvectors = eigenvectors(:, id);
    sigma_d_evals_sorted = sigma_d_evals(id);
    c2_oo_evals_sorted = c2_oo_evals(id);
    inv_evals_sorted = inv_evals(id);
    c2_o_evals_sorted = c2_o_evals(id);
    
    % Analyze for degeneracies
    energy = energies(1);
    degeneracy = 1;
    index = 1;
    trace_sigma_d = sigma_d_evals_sorted(1);
    trace_c2_oo = c2_oo_evals_sorted(1);
    trace_inv = inv_evals_sorted(1);
    trace_c2_o = c2_o_evals_sorted(1);
    
    for i = 2:(total_sites)
        if (abs(energies(i) - energy) < tolerance) % Next energy is degenerate
            degeneracy = degeneracy + 1;
            trace_sigma_d = trace_sigma_d + sigma_d_evals_sorted(i);
            trace_c2_oo = trace_c2_oo + c2_oo_evals_sorted(i);
            trace_inv = trace_inv + inv_evals_sorted(i);
            trace_c2_o = trace_c2_o + c2_o_evals_sorted(i);
        else % Next energy is new
            % Record stats for previous energy level
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_c2_oo;
            energy_levels(index, 4) = trace_c2_o;
            energy_levels(index, 5) = trace_inv;
            energy_levels(index, 6) = trace_sigma_d;
            energy_levels(index, 7) = i - 1;
    
            % Reset variables for new energy level
            energy = energies(i);
            degeneracy = 1;
            index = index + 1;
            trace_c2_oo = c2_oo_evals_sorted(i);
            trace_c2_o = c2_o_evals_sorted(i);
            trace_inv = inv_evals_sorted(i);
            trace_sigma_d = sigma_d_evals_sorted(i);
        end
    
        % Record energy level if we reach the end
        if (i == total_sites)
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_c2_oo;
            energy_levels(index, 4) = trace_c2_o;
            energy_levels(index, 5) = trace_inv;
            energy_levels(index, 6) = trace_sigma_d;
            energy_levels(index, 7) = i;
            
        end
    end
    
    energy_levels = energy_levels(1:index, :);
    
    %% Separate the representation classes
    % Allocate space for each irrep
    elevels_a1g = zeros(total_sites, 7);
    elevels_a2g = zeros(total_sites, 7);
    elevels_b1g = zeros(total_sites, 7);
    elevels_b2g = zeros(total_sites, 7);
    elevels_eg = zeros(total_sites, 7);
    elevels_a1u = zeros(total_sites, 7);
    elevels_a2u = zeros(total_sites, 7);
    elevels_b1u = zeros(total_sites, 7);
    elevels_b2u = zeros(total_sites, 7);
    elevels_eu = zeros(total_sites, 7);

    % Allocate space for solvable elevels
    elevels_a1g_solv = zeros(total_sites, 7);
    elevels_a1u_solv = zeros(total_sites, 7);

    if (l1_even == false && l2_even == true)
        elevels_b2g_solv = zeros(total_sites, 7);
        elevels_b2u_solv = zeros(total_sites, 7);
    elseif (l1_even == true && l2_even == false)
        elevels_a2g_solv = zeros(total_sites, 7);
        elevels_a2u_solv = zeros(total_sites, 7);
    elseif (l1_even == false && l2_even == false)
        elevels_b1g_solv = zeros(total_sites, 7);
        elevels_b1u_solv = zeros(total_sites, 7);
    end

    % Set indices
    index_a1g = 1;
    index_a2g = 1;
    index_b1g = 1;
    index_b2g = 1;
    index_eg = 1;
    index_a1u = 1;
    index_a2u = 1;
    index_b1u = 1;
    index_b2u = 1;
    index_eu = 1;

    index_a1g_solv = 1;
    index_a1u_solv = 1;

    if (l1_even == false && l2_even == true)
        index_b2g_solv = 1;
        index_b2u_solv = 1;
    elseif (l1_even == true && l2_even == false)
        index_a2g_solv = 1;
        index_a2u_solv = 1;
    elseif (l1_even == false && l2_even == false)
        index_b1g_solv = 1;
        index_b1u_solv = 1;
    end

    % Round the energy_levels characters
    energy_levels_rounded = zeros(size(energy_levels, 1), 7);
    energy_levels_rounded(:, 2:6) = round(energy_levels(:, 2:6));
    energy_levels_rounded(:, 1) = energy_levels(:, 1);
    energy_levels_rounded(:, 7) = energy_levels(:, 7);
    
    % Fill in elevel info in each irrep
    for i = 1:size(energy_levels_rounded, 1)
        trace_e = energy_levels_rounded(i, 2);
        trace_c2_oo = energy_levels_rounded(i, 3);
        trace_c2_o = energy_levels_rounded(i, 4);
        trace_inv = energy_levels_rounded(i, 5);
        trace_sigma_d = energy_levels_rounded(i, 6);
    
        traces = [trace_e, trace_c2_oo, trace_c2_o, trace_inv, trace_sigma_d];
    
        if (isequal(traces, [1, 1, 1, 1, 1])) % A1g
            elevels_a1g(index_a1g, :) = energy_levels_rounded(i, :);
            index_a1g = index_a1g + 1;
        elseif (isequal(traces, [2, 2, 2, 0, 0])) % Accidental degeneracy A1g + A1u
            elevels_a1g_solv(index_a1g_solv, :) = energy_levels_rounded(i, :);
            index_a1g_solv = index_a1g_solv + 1;
            elevels_a1u_solv(index_a1u_solv, :) = energy_levels_rounded(i, :); 
            index_a1u_solv = index_a1u_solv + 1;
        elseif (isequal(traces, [2, -2, -2, 0, 0])) % Accidental degeneracy A2g + A2u
            elevels_a2g_solv(index_a2g_solv, :) = energy_levels_rounded(i, :);
            index_a2g_solv = index_a2g_solv + 1;
            elevels_a2u_solv(index_a2u_solv, :) = energy_levels_rounded(i, :); 
            index_a2u_solv = index_a2u_solv + 1;
        elseif (isequal(traces, [2, 2, -2, 0, 0])) % Accidental degeneracy B1g + B1u
            elevels_b1g_solv(index_b1g_solv, :) = energy_levels_rounded(i, :);
            index_b1g_solv = index_b1g_solv + 1;
            elevels_b1u_solv(index_b1u_solv, :) = energy_levels_rounded(i, :); 
            index_b1u_solv = index_b1u_solv + 1;
        elseif (isequal(traces, [2, -2, 2, 0, 0])) % Accidental degeneracy B2g + B2u
            elevels_b2g_solv(index_b2g_solv, :) = energy_levels_rounded(i, :);
            index_b2g_solv = index_b2g_solv + 1;
            elevels_b2u_solv(index_b2u_solv, :) = energy_levels_rounded(i, :);
            index_b2u_solv = index_b2u_solv + 1;
        elseif (isequal(traces, [1, -1, -1, 1, -1])) % A2g
            elevels_a2g(index_a2g, :) = energy_levels_rounded(i, :);
            index_a2g = index_a2g + 1;
        elseif (isequal(traces, [1, 1, -1, 1, -1])) % B1g
            elevels_b1g(index_b1g, :) = energy_levels_rounded(i, :);
            index_b1g = index_b1g + 1;
        elseif (isequal(traces, [1, -1, 1, 1, 1])) % B2g
            elevels_b2g(index_b2g, :) = energy_levels_rounded(i, :);
            index_b2g = index_b2g + 1;
        elseif (isequal(traces, [2, 0, 0, 2, 0])) % Eg
            elevels_eg(index_eg, :) = energy_levels_rounded(i, :);
            index_eg = index_eg + 1;
        elseif (isequal(traces, [1, 1, 1, -1, -1])) % A1u
            elevels_a1u(index_a1u, :) = energy_levels_rounded(i, :);
            index_a1u = index_a1u + 1;
        elseif (isequal(traces, [1, -1, -1, -1, 1])) % A2u
            elevels_a2u(index_a2u, :) = energy_levels_rounded(i, :);
            index_a2u = index_a2u + 1;
        elseif (isequal(traces, [1, 1, -1, -1, 1])) % B1u
            elevels_b1u(index_b1u, :) = energy_levels_rounded(i, :);
            index_b1u = index_b1u + 1;
        elseif (isequal(traces, [1, -1, 1, -1, -1])) % B2u
            elevels_b2u(index_b2u, :) = energy_levels_rounded(i, :);
            index_b2u = index_b2u + 1;
        elseif (isequal(traces, [2, 0, 0, -2, 0])) % Eu
            elevels_eu(index_eu, :) = energy_levels_rounded(i, :);
            index_eu = index_eu + 1;
        else 
            fprintf([num2str(energy_levels_rounded(i, 1)) ' ' num2str(energy_levels_rounded(i, 2)) ' ' ...
                num2str(energy_levels_rounded(i, 3)) ' ' num2str(energy_levels_rounded(i, 4)) ' ' ...
                num2str(energy_levels_rounded(i, 5)) ' ' num2str(energy_levels_rounded(i, 6)) ' ' ...
                num2str(energy_levels_rounded(i, 7)) '\n'])
        end
    end

    %% Calculate energies of non-degen solvable wavefunctions
    % Calculate solvable wavefunctions energies from A1g (n1, n1) with n1 even
    a1g_exp_energies_even = zeros(N + 1, 1);
    for i = 0:2:(N - 1)
        n1 = i;
        a1g_exp_energies_even(i + 1) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
        a1g_exp_energies_even(i + 2) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * 0 / N));
    end

    if (l1_even == false && l2_even) % Calculate solvable wavefunctions energies from B2g (n1, n1) with n1 odd (l1 odd, l2 even)
        b2g_exp_energies_odd = zeros(ceil(N / 2), 1);
        for i = 1:ceil(N / 2)
            n1 = (i - 1) * 2 + 1;
            b2g_exp_energies_odd(i) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
        end
    elseif (l1_even && l2_even == false) % Calculate solvable wavefunctions energies from A2u (n1, n1) with n1 odd (l1 even, l2 odd)
        a2u_exp_energies_odd = zeros(ceil(N / 2), 1);
        for i = 1:ceil(N / 2)
            n1 = (i - 1) * 2 + 1;
            a2u_exp_energies_odd(i) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
        end
    elseif (l1_even == false && l2_even == false) % Calculate solvable wavefunctions energies from A2u (n1, n1) with n1 odd (l1 even, l2 odd)
        b1u_exp_energies_odd = zeros(ceil(N / 2), 1);
        for i = 1:ceil(N / 2)
            n1 = (i - 1) * 2 + 1;
            b1u_exp_energies_odd(i) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
        end
    end

    % Calculate edge number wave energies (N, n1)
    edge_wave_energies = zeros(ceil(N / 2), 1);
    for i = 1:ceil(N / 2)
        n1 = 2 * i - 1;
        edge_wave_energies(i) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * N / N));
    end

    %% Remove non-degen and edge solvable wavefunctions

    % Remove from A1g (even n1)
    for i = 1:(index_a1g - 1)
        if (i < index_a1g)
            for j = 1:(N + 1)
                if (abs(a1g_exp_energies_even(j) - elevels_a1g(i, 1)) < tolerance)
                    % Add to elevel_a1g_solv
                    elevels_a1g_solv(index_a1g_solv, :) = elevels_a1g(i, :);
                    index_a1g_solv = index_a1g_solv + 1;
    
                    % Delete from elevels_a1g
                    elevels_a1g(i, :) = [];
                    index_a1g = index_a1g - 1;
                end
            end 
        end
    end

    if (l1_even == false && l2_even) % Remove from B2g (odd n1)
        for i = 1:(index_b2g - 1)
            if (i < index_b2g)
                for j = 1:ceil(N / 2)
                    if (abs(b2g_exp_energies_odd(j) - elevels_b2g(i, 1)) < tolerance)
                        % Add to elevel_b2g_solv
                        elevels_b2g_solv(index_b2g_solv, :) = elevels_b2g(i, :);
                        index_b2g_solv = index_b2g_solv + 1;
    
                        % Delete from elevels_b2g
                        elevels_b2g(i, :) = [];
                        index_b2g = index_b2g - 1;
                    end

                    if ((i < index_b2g) && (abs(edge_wave_energies(j) - elevels_b2g(i, 1)) < tolerance))
                        % Add to elevel_b2g_solv
                        elevels_b2g_solv(index_b2g_solv, :) = elevels_b2g(i, :);
                        index_b2g_solv = index_b2g_solv + 1;
    
                        % Delete from elevels_b2g
                        elevels_b2g(i, :) = [];
                        index_b2g = index_b2g - 1;
                    end
                end
            end
        end
        for i = 1:(index_b2u - 1)
            if (i < index_b2u)
                for j = 1:ceil(N / 2)
                    if ((i < index_b2u) && (abs(edge_wave_energies(j) - elevels_b2u(i, 1)) < tolerance))
                        % Add to elevel_b2u_solv
                        elevels_b2u_solv(index_b2u_solv, :) = elevels_b2u(i, :);
                        index_b2u_solv = index_b2u_solv + 1;
    
                        % Delete from elevels_b2u
                        elevels_b2u(i, :) = [];
                        index_b2u = index_b2u - 1;
                    end
                end
            end
        end      
    elseif (l1_even && l2_even == false) % Remove from A2u (odd n1)
        for i = 1:(index_a2u - 1)
            if (i < index_a2u)
                for j = 1:ceil(N / 2)
                    if (abs(a2u_exp_energies_odd(j) - elevels_a2u(i, 1)) < tolerance)
                        % Add to elevel_a2u_solv
                        elevels_a2u_solv(index_a2u_solv, :) = elevels_a2u(i, :);
                        index_a2u_solv = index_a2u_solv + 1;
    
                        % Delete from elevels_a2u
                        elevels_a2u(i, :) = [];
                        index_a2u = index_a2u - 1;
                    end

                    if ((i < index_a2u) && (abs(edge_wave_energies(j) - elevels_a2u(i, 1)) < tolerance))
                        % Add to elevel_a2u_solv
                        elevels_a2u_solv(index_a2u_solv, :) = elevels_a2u(i, :);
                        index_a2u_solv = index_a2u_solv + 1;
    
                        % Delete from elevels_a2u
                        elevels_a2u(i, :) = [];
                        index_a2u = index_a2u - 1;
                    end
                end
            end
        end
        for i = 1:(index_a2g - 1)
            if (i < index_a2g)
                for j = 1:ceil(N / 2)
                    if ((i < index_a2g) && (abs(edge_wave_energies(j) - elevels_a2g(i, 1)) < tolerance))
                        % Add to elevel_a2g_solv
                        elevels_a2g_solv(index_a2g_solv, :) = elevels_a2g(i, :);
                        index_a2g_solv = index_a2g_solv + 1;
    
                        % Delete from elevels_a2g
                        elevels_a2g(i, :) = [];
                        index_a2g = index_a2g - 1;
                    end
                end
            end
        end 
    elseif (l1_even == false && l2_even == false) % Remove from B1u (odd n1)
        for i = 1:(index_b1u - 1)
            if (i < index_b1u)
                for j = 1:ceil(N / 2)
                    if (abs(b1u_exp_energies_odd(j) - elevels_b1u(i, 1)) < tolerance)
                        % Add to elevel_b1u_solv
                        elevels_b1u_solv(index_b1u_solv, :) = elevels_b1u(i, :);
                        index_b1u_solv = index_b1u_solv + 1;
    
                        % Delete from elevels_b1u
                        elevels_b1u(i, :) = [];
                        index_b1u = index_b1u - 1;
                    end

                    if ((i < index_b1u) && (abs(edge_wave_energies(j) - elevels_b1u(i, 1)) < tolerance))
                        % Add to elevel_b1u_solv
                        elevels_b1u_solv(index_b1u_solv, :) = elevels_b1u(i, :);
                        index_b1u_solv = index_b1u_solv + 1;
    
                        % Delete from elevels_b1u
                        elevels_b1u(i, :) = [];
                        index_b1u = index_b1u - 1;
                    end
                end
            end
        end
        for i = 1:(index_b1g - 1)
            if (i < index_b1g)
                for j = 1:ceil(N / 2)
                    if ((i < index_b1g) && (abs(edge_wave_energies(j) - elevels_b1g(i, 1)) < tolerance))
                        % Add to elevel_b1g_solv
                        elevels_b1g_solv(index_b1g_solv, :) = elevels_b1g(i, :);
                        index_b1g_solv = index_b1g_solv + 1;
    
                        % Delete from elevels_b1g
                        elevels_b1g(i, :) = [];
                        index_b1g = index_b1g - 1;
                    end
                end
            end
        end 
    end
 
    %% Clean up
    % Remove extra rows of zeros
    if (index_a1g > 1)
        elevels_a1g = elevels_a1g(1:(index_a1g - 1), :);
    end 
    
    if (index_a2g > 1)
        elevels_a2g = elevels_a2g(1:(index_a2g - 1), :);
    end
    
    if (index_b1g > 1)
        elevels_b1g = elevels_b1g(1:(index_b1g - 1), :);
    end
    
    if (index_b2g > 1)
        elevels_b2g = elevels_b2g(1:(index_b2g - 1), :);
    end
    
    if (index_eg > 1)
        elevels_eg = elevels_eg(1:(index_eg - 1), :);
    end
    
    if (index_a1u > 1)
        elevels_a1u = elevels_a1u(1:(index_a1u - 1), :);
    end 
    
    if (index_a2u > 1)
        elevels_a2u = elevels_a2u(1:(index_a2u - 1), :);
    end
    
    if (index_b1u > 1)
        elevels_b1u = elevels_b1u(1:(index_b1u - 1), :);
    end
    
    if (index_b2u > 1)
        elevels_b2u = elevels_b2u(1:(index_b2u - 1), :);
    end
    
    if (index_eu > 1)
        elevels_eu = elevels_eu(1:(index_eu - 1), :);
    end

    if (index_a1g_solv > 1)
        elevels_a1g_solv = elevels_a1g_solv(1:(index_a1g_solv - 1), :);
    end 

    if (index_a1u_solv > 1)
        elevels_a1u_solv = elevels_a1u_solv(1:(index_a1u_solv - 1), :);
    end 
    
    if (l1_even == false && l2_even)
        if (index_b2g_solv > 1)
            elevels_b2g_solv = elevels_b2g_solv(1:(index_b2g_solv - 1), :);
        end

        if (index_b2u_solv > 1)
            elevels_b2u_solv = elevels_b2u_solv(1:(index_b2u_solv - 1), :);
        end
    elseif (l1_even && l2_even == false)
        if (index_a2g_solv > 1)
            elevels_a2g_solv = elevels_a2g_solv(1:(index_a2g_solv - 1), :);
        end

        if (index_a2u_solv > 1)
            elevels_a2u_solv = elevels_a2u_solv(1:(index_a2u_solv - 1), :);
        end
    elseif (l1_even == false && l2_even == false)
        if (index_b1g_solv > 1)
            elevels_b1g_solv = elevels_b1g_solv(1:(index_b1g_solv - 1), :);
        end

        if (index_b1u_solv > 1)
            elevels_b1u_solv = elevels_b1u_solv(1:(index_b1u_solv - 1), :);
        end
    end

    % Resort elevels_solv
    elevels_a1g_solv_sorted = sortrows(elevels_a1g_solv);
    elevels_a1u_solv_sorted = sortrows(elevels_a1u_solv);

    if (l1_even == false && l2_even)  
        elevels_b2g_solv_sorted = sortrows(elevels_b2g_solv);
        elevels_b2u_solv_sorted = sortrows(elevels_b2u_solv);
    elseif (l1_even && l2_even == false)  
        elevels_a2u_solv_sorted = sortrows(elevels_a2u_solv);
        elevels_a2g_solv_sorted = sortrows(elevels_a2g_solv);
    elseif (l1_even == false && l2_even == false)  
        elevels_b1u_solv_sorted = sortrows(elevels_b1u_solv);
        elevels_b1g_solv_sorted = sortrows(elevels_b1g_solv);
    end

    %% Find energy level spacings STOPPED HERE
    % Make vector of energy level spacings
    spacings_a1g = zeros(size(elevels_a1g, 1) - 1, 1);
    spacings_a2g = zeros(size(elevels_a2g, 1) - 1, 1);
    spacings_b1g = zeros(size(elevels_b1g, 1) - 1, 1);
    spacings_b2g = zeros(size(elevels_b2g, 1) - 1, 1);
    spacings_eg = zeros(size(elevels_eg, 1) - 1, 1);
    spacings_a1u = zeros(size(elevels_a1u, 1) - 1, 1);
    spacings_a2u = zeros(size(elevels_a2u, 1) - 1, 1);
    spacings_b1u = zeros(size(elevels_b1u, 1) - 1, 1);
    spacings_b2u = zeros(size(elevels_b2u, 1) - 1, 1);
    spacings_eu = zeros(size(elevels_eu, 1) - 1, 1);

    spacings_a1g_solv = zeros(size(elevels_a1g_solv, 1) - 1, 1);
    spacings_a1u_solv = zeros(size(elevels_a1u_solv, 1) - 1, 1);

    if (l1_even == false && l2_even)
        spacings_b2g_solv = zeros(size(elevels_b2g_solv, 1) - 1, 1);
        spacings_b2u_solv = zeros(size(elevels_b2u_solv, 1) - 1, 1);
    elseif (l1_even && l2_even == false)
        spacings_a2g_solv = zeros(size(elevels_a2g_solv, 1) - 1, 1);
        spacings_a2u_solv = zeros(size(elevels_a2u_solv, 1) - 1, 1);
    elseif (l1_even == false && l2_even == false)
        spacings_b1g_solv = zeros(size(elevels_b1g_solv, 1) - 1, 1);
        spacings_b1u_solv = zeros(size(elevels_b1u_solv, 1) - 1, 1);
    end
    
    for i = 1:(size(elevels_a1g, 1) - 1)
        spacings_a1g(i) = abs(elevels_a1g(i, 1) - elevels_a1g(i+1, 1));
    end
    
    for i = 1:(size(elevels_a2g, 1) - 1)
        spacings_a2g(i) = abs(elevels_a2g(i, 1) - elevels_a2g(i+1, 1));
    end
    
    for i = 1:(size(elevels_b1g, 1) - 1)
        spacings_b1g(i) = abs(elevels_b1g(i, 1) - elevels_b1g(i+1, 1));
    end
    
    for i = 1:(size(elevels_b2g, 1) - 1)
        spacings_b2g(i) = abs(elevels_b2g(i, 1) - elevels_b2g(i+1, 1));
    end
    
    for i = 1:(size(elevels_eg, 1) - 1)
        spacings_eg(i) = abs(elevels_eg(i, 1) - elevels_eg(i+1, 1));
    end
    
    for i = 1:(size(elevels_a1u, 1) - 1)
        spacings_a1u(i) = abs(elevels_a1u(i, 1) - elevels_a1u(i+1, 1));
    end
    
    for i = 1:(size(elevels_a2u, 1) - 1)
        spacings_a2u(i) = abs(elevels_a2u(i, 1) - elevels_a2u(i+1, 1));
    end
    
    for i = 1:(size(elevels_b1u, 1) - 1)
        spacings_b1u(i) = abs(elevels_b1u(i, 1) - elevels_b1u(i+1, 1));
    end
    
    for i = 1:(size(elevels_b2u, 1) - 1)
        spacings_b2u(i) = abs(elevels_b2u(i, 1) - elevels_b2u(i+1, 1));
    end
    
    for i = 1:(size(elevels_eu, 1) - 1)
        spacings_eu(i) = abs(elevels_eu(i, 1) - elevels_eu(i+1, 1));
    end

    for i = 1:(size(elevels_a1g_solv_sorted, 1) - 1)
        spacings_a1g_solv(i) = abs(elevels_a1g_solv_sorted(i, 1) - elevels_a1g_solv_sorted(i + 1, 1));
    end

    for i = 1:(size(elevels_a1u_solv_sorted, 1) - 1)
        spacings_a1u_solv(i) = abs(elevels_a1u_solv_sorted(i, 1) - elevels_a1u_solv_sorted(i+1, 1));
    end

    if (l1_even == false && l2_even)
        for i = 1:(size(elevels_b2g_solv_sorted, 1) - 1)
            spacings_b2g_solv(i) = abs(elevels_b2g_solv_sorted(i, 1) - elevels_b2g_solv_sorted(i+1, 1));
        end
        
        for i = 1:(size(elevels_b2u_solv_sorted, 1) - 1)
            spacings_b2u_solv(i) = abs(elevels_b2u_solv_sorted(i, 1) - elevels_b2u_solv_sorted(i+1, 1));
        end
    elseif (l1_even && l2_even == false)
        for i = 1:(size(elevels_a2u_solv_sorted, 1) - 1)
            spacings_a2u_solv(i) = abs(elevels_a2u_solv_sorted(i, 1) - elevels_a2u_solv_sorted(i+1, 1));
        end
        
        for i = 1:(size(elevels_a2g_solv_sorted, 1) - 1)
            spacings_a2g_solv(i) = abs(elevels_a2g_solv_sorted(i, 1) - elevels_a2g_solv_sorted(i+1, 1));
        end
    elseif (l1_even == false && l2_even == false)
        for i = 1:(size(elevels_b1u_solv_sorted, 1) - 1)
            spacings_b1u_solv(i) = abs(elevels_b1u_solv_sorted(i, 1) - elevels_b1u_solv_sorted(i+1, 1));
        end
        
        for i = 1:(size(elevels_b1g_solv_sorted, 1) - 1)
            spacings_b1g_solv(i) = abs(elevels_b1g_solv_sorted(i, 1) - elevels_b1g_solv_sorted(i+1, 1));
        end
    end
    
    % Compute level spacing ratios
    r_a1g = LSR(spacings_a1g);
    r_a2g = LSR(spacings_a2g);
    r_b1g = LSR(spacings_b1g);
    r_b2g = LSR(spacings_b2g);
    r_eg = LSR(spacings_eg);
    r_a1u = LSR(spacings_a1u);
    r_a2u = LSR(spacings_a2u);
    r_b1u = LSR(spacings_b1u);
    r_b2u = LSR(spacings_b2u);
    r_eu = LSR(spacings_eu);

    r_a1g_solv = LSR(spacings_a1g_solv);
    r_a1u_solv = LSR(spacings_a1u_solv);
    
    if (l1_even == false && l2_even)
        r_b2g_solv = LSR(spacings_b2g_solv);
        r_b2u_solv = LSR(spacings_b2u_solv);
    elseif (l1_even && l2_even == false)
        r_a2g_solv = LSR(spacings_a2g_solv);
        r_a2u_solv = LSR(spacings_a2u_solv);
    elseif (l1_even == false && l2_even == false)
        r_b1g_solv = LSR(spacings_b1g_solv);
        r_b1u_solv = LSR(spacings_b1u_solv);
    end

    % Fill in r array values
    if (l1_even == false && l2_even)
        r_array(index_r, :) = [N, r_a1g, r_a2g, r_b1g, r_b2g, r_eg, r_a1u, r_a2u, r_b1u, r_b2u, r_eu, r_a1g_solv, r_a1u_solv, r_b2g_solv, r_b2u_solv];
    elseif (l1_even && l2_even == false)
        r_array(index_r, :) = [N, r_a1g, r_a2g, r_b1g, r_b2g, r_eg, r_a1u, r_a2u, r_b1u, r_b2u, r_eu, r_a1g_solv, r_a1u_solv, r_a2g_solv, r_a2u_solv];
    elseif (l1_even == false && l2_even == false)
        r_array(index_r, :) = [N, r_a1g, r_a2g, r_b1g, r_b2g, r_eg, r_a1u, r_a2u, r_b1u, r_b2u, r_eu, r_a1g_solv, r_a1u_solv, r_b1g_solv, r_b1u_solv];
    end

    % Fill in size of irreps
    if (l1_even == false && l2_even)
        size_array(index_size, :) = [N, size(elevels_a1g, 1), size(elevels_a2g, 1), size(elevels_b1g, 1), ...
            size(elevels_b2g, 1), size(elevels_eg, 1), size(elevels_a1u, 1), size(elevels_a2u, 1), ...
            size(elevels_b1u, 1), size(elevels_b2u, 1), size(elevels_eu, 1), size(elevels_a1g_solv, 1), ...
            size(elevels_a1u_solv, 1), size(elevels_b2g_solv, 1), size(elevels_b2u_solv, 1)];
    elseif (l1_even && l2_even == false)
        size_array(index_size, :) = [N, size(elevels_a1g, 1), size(elevels_a2g, 1), size(elevels_b1g, 1), ...
            size(elevels_b2g, 1), size(elevels_eg, 1), size(elevels_a1u, 1), size(elevels_a2u, 1), ...
            size(elevels_b1u, 1), size(elevels_b2u, 1), size(elevels_eu, 1), size(elevels_a1g_solv, 1), ...
            size(elevels_a1u_solv, 1), size(elevels_a2g_solv, 1), size(elevels_a2u_solv, 1)];
    elseif (l1_even == false && l2_even == false)
        size_array(index_size, :) = [N, size(elevels_a1g, 1), size(elevels_a2g, 1), size(elevels_b1g, 1), ...
            size(elevels_b2g, 1), size(elevels_eg, 1), size(elevels_a1u, 1), size(elevels_a2u, 1), ...
            size(elevels_b1u, 1), size(elevels_b2u, 1), size(elevels_eu, 1), size(elevels_a1g_solv, 1), ...
            size(elevels_a1u_solv, 1), size(elevels_b1g_solv, 1), size(elevels_b1u_solv, 1)];
    end

    % Find proportion of solvable states
    if (l1_even == false && l2_even) % l1 is odd
        solvable_prop(index_size, :) = [size(elevels_a1g_solv, 1), size(elevels_a1u_solv, 1), size(elevels_b2g_solv, 1), size(elevels_b2u_solv, 1)] / total_sites;
    elseif (l1_even && l2_even == false) % l1 is even
        solvable_prop(index_size, :) = [size(elevels_a1g_solv, 1), size(elevels_a1u_solv, 1), size(elevels_a2g_solv, 1), size(elevels_a2u_solv, 1)] / total_sites;
    elseif (l1_even == false && l2_even == false)
        solvable_prop(index_size, :) = [size(elevels_a1g_solv, 1), size(elevels_a1u_solv, 1), size(elevels_b1g_solv, 1), size(elevels_b1u_solv, 1)] / total_sites;
    end

    %% Plot the level spacings histogram and show r values
    % Plot histograms
    figure(index_n)
    tiledlayout(2,5,'TileSpacing', 'tight','Padding','Tight')
    bin_factor = 5;
    
    nexttile
    histogram(spacings_a1g, ceil(size(elevels_a1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A1g')
    subtitle(['r = ' num2str(r_a1g) '; total = ' num2str(size_array(index_size, 2))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_a2g, ceil(size(elevels_a2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A2g')
    subtitle(['r = ' num2str(r_a2g) '; total = ' num2str(size_array(index_size, 3))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b1g, ceil(size(elevels_b1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B1g')
    subtitle(['r = ' num2str(r_b1g) '; total = ' num2str(size_array(index_size, 4))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b2g, ceil(size(elevels_b2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B2g')
    subtitle(['r = ' num2str(r_b2g) '; total = ' num2str(size_array(index_size, 5))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_eg, ceil(size(elevels_eg, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Eg')
    subtitle(['r = ' num2str(r_eg) '; total = ' num2str(size_array(index_size, 6))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_a1u, ceil(size(elevels_a1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A1u')
    subtitle(['r = ' num2str(r_a1u) '; total = ' num2str(size_array(index_size, 7))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_a2u, ceil(size(elevels_a2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A2u')
    subtitle(['r = ' num2str(r_a2u) '; total = ' num2str(size_array(index_size, 8))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b1u, ceil(size(elevels_b1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B1u')
    subtitle(['r = ' num2str(r_b1u) '; total = ' num2str(size_array(index_size, 9))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b2u, ceil(size(elevels_b2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B2u')
    subtitle(['r = ' num2str(r_b2u) '; total = ' num2str(size_array(index_size, 10))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_eu, ceil(size(elevels_eu, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Eu')
    subtitle(['r = ' num2str(r_eu) '; total = ' num2str(size_array(index_size, 11))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    set(figure(index_n),'position',[0,100,2000,400])
    saveas(gcf, [folderName '/N=' num2str(N) '_elevel_hist_unsolv.jpeg']);

    % Plot histogram of solvable wavefunctions
    figure(index_n + 1)
        tiledlayout(2,2, 'TileSpacing', 'tight','Padding','Tight')
        bin_factor = 5;

        nexttile
        histogram(spacings_a1g_solv, ceil(size(elevels_a1g_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('A1g Solvable')
        subtitle(['r = ' num2str(r_a1g_solv) '; total = ' num2str(size(elevels_a1g_solv, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    
        nexttile
        histogram(spacings_a1u_solv, ceil(size(elevels_a1u_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('A1u Solvable')
        subtitle(['r = ' num2str(r_a1u_solv) '; total = ' num2str(size(elevels_a1u_solv, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')

    if (l1_even == false && l2_even)
        nexttile
        histogram(spacings_b2g_solv, ceil(size(elevels_b2g_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('B2g Solvable')
        subtitle(['r = ' num2str(r_b2g_solv) '; total = ' num2str(size(elevels_b2g_solv, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    
        nexttile
        histogram(spacings_b2u_solv, ceil(size(elevels_b2u_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('B2u Solvable')
        subtitle(['r = ' num2str(r_b2u_solv) '; total = ' num2str(size(elevels_b2u_solv, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    elseif (l1_even && l2_even == false)
        nexttile
        histogram(spacings_a2g_solv, ceil(size(elevels_a2g_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('A2g Solvable')
        subtitle(['r = ' num2str(r_a2g_solv) '; total = ' num2str(size(elevels_a2g_solv, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    
        nexttile
        histogram(spacings_a2u_solv, ceil(size(elevels_a2u_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('A2u Solvable')
        subtitle(['r = ' num2str(r_a2u_solv) '; total = ' num2str(size(elevels_a2u_solv, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    elseif (l1_even == false && l2_even == false)
        nexttile
        histogram(spacings_b1g_solv, ceil(size(elevels_b1g_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('B1g Solvable')
        subtitle(['r = ' num2str(r_b1g_solv) '; total = ' num2str(size(elevels_b1g_solv, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    
        nexttile
        histogram(spacings_b1u_solv, ceil(size(elevels_b1u_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('B1u Solvable')
        subtitle(['r = ' num2str(r_b1u_solv) '; total = ' num2str(size(elevels_b1u_solv, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    end
    set(figure(index_n + 1),'position',[0, 100, 700, 400])
    saveas(gcf, [folderName '/N=' num2str(N) '_elevel_hist_solv.jpeg']);


    % Plot histogram after removal of solvable wavefunctions
    figure(index_n + 2)
        tiledlayout(2,2, 'TileSpacing', 'tight','Padding','Tight')
        bin_factor = 5;

        nexttile
        histogram(spacings_a1g, ceil(size(elevels_a1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('A1g')
        subtitle(['r = ' num2str(r_a1g) '; total = ' num2str(size(elevels_a1g, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    
        nexttile
        histogram(spacings_a1u, ceil(size(elevels_a1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('A1u')
        subtitle(['r = ' num2str(r_a1u) '; total = ' num2str(size(elevels_a1u, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')

    if (l1_even == false && l2_even)
        nexttile
        histogram(spacings_b2g, ceil(size(elevels_b2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('B2g')
        subtitle(['r = ' num2str(r_b2g) '; total = ' num2str(size(elevels_b2g, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    
        nexttile
        histogram(spacings_b2u, ceil(size(elevels_b2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('B2u')
        subtitle(['r = ' num2str(r_b2u) '; total = ' num2str(size(elevels_b2u, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    elseif (l1_even && l2_even == false)
        nexttile
        histogram(spacings_a2g, ceil(size(elevels_a2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('A2g')
        subtitle(['r = ' num2str(r_a2g) '; total = ' num2str(size(elevels_a2g, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    
        nexttile
        histogram(spacings_a2u, ceil(size(elevels_a2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('A2u')
        subtitle(['r = ' num2str(r_a2u) '; total = ' num2str(size(elevels_a2u, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    elseif (l1_even == false && l2_even == false)
        nexttile
        histogram(spacings_b1g, ceil(size(elevels_b1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('B1g')
        subtitle(['r = ' num2str(r_b1g) '; total = ' num2str(size(elevels_b1g, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    
        nexttile
        histogram(spacings_b1u, ceil(size(elevels_b1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('B1u')
        subtitle(['r = ' num2str(r_b1u) '; total = ' num2str(size(elevels_b1u, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    end
    set(figure(index_n + 2),'position',[0, 100, 700, 400])
    saveas(gcf, [folderName '/N=' num2str(N) '_elevel_hist_unsolv_removed.jpeg']);

    % Increment indices
    index_n = index_n + 3;
    index_size = index_size + 1;
    index_r = index_r + 1;
    
end

%% Plot r as a function of size in each irrep
figure(index_n)
tiledlayout(2,5,'TileSpacing', 'tight','Padding','Tight')
ylow_unsolv = 0.45;
yhigh_unsolv = 0.7;
ylow_solv = 0.35;
yhigh_solv = 0.6;

nexttile
plot(N_nums.', r_array(:, 2).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A1g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 3).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A2g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 4).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B1g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 5).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B2g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 6).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Eg')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 7).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A1u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 8).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A2u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 9).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B1u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 10).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('B2u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 11).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Eu')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

set(figure(index_n),'position',[0,100,3000,400])
saveas(gcf, [folderName '/LSR_plot.jpeg']);

% Plot r values for solvable wavefunctions
figure(index_n + 1)
    tiledlayout(2, 2, 'TileSpacing', 'tight','Padding','Tight')

    nexttile
    plot(N_nums.', r_array(:, 12).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    title('A1g')
    xlabel('N')
    ylabel('r')
    hold on 
    yline(0.39, '--', 'r = 0.39')
    ylim([ylow_solv, yhigh_solv])

    nexttile
    plot(N_nums, r_array(:, 13).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    title('A1u')
    xlabel('N')
    ylabel('r')
    hold on 
    yline(0.39, '--', 'r = 0.39')
    ylim([ylow_solv, yhigh_solv])

if (l1_even == false && l2_even)
    nexttile
    plot(N_nums, r_array(:, 14).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    title('B2g')
    xlabel('N')
    ylabel('r')
    hold on 
    yline(0.39, '--', 'r = 0.39')
    ylim([ylow_solv, yhigh_solv])
    
    nexttile
    plot(N_nums, r_array(:, 15).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    title('B2u')
    xlabel('N')
    ylabel('r')
    hold on 
    yline(0.39, '--', 'r = 0.39')
    ylim([ylow_solv, yhigh_solv])
elseif (l1_even && l2_even == false)
    nexttile
    plot(N_nums, r_array(:, 14).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    title('A2g')
    xlabel('N')
    ylabel('r')
    hold on 
    yline(0.39, '--', 'r = 0.39')
    ylim([ylow_solv, yhigh_solv])
    
    nexttile
    plot(N_nums, r_array(:, 15).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    title('A2u')
    xlabel('N')
    ylabel('r')
    hold on 
    yline(0.39, '--', 'r = 0.39')
    ylim([ylow_solv, yhigh_solv])
elseif (l1_even == false && l2_even == false)
    nexttile
    plot(N_nums, r_array(:, 14).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    title('B1g')
    xlabel('N')
    ylabel('r')
    hold on 
    yline(0.39, '--', 'r = 0.39')
    ylim([ylow_solv, yhigh_solv])
    
    nexttile
    plot(N_nums, r_array(:, 15).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    title('B1u')
    xlabel('N')
    ylabel('r')
    hold on 
    yline(0.39, '--', 'r = 0.39')
    ylim([ylow_solv, yhigh_solv])
end

set(figure(index_n + 1),'position',[0,100,1000,600])
saveas(gcf, [folderName '/LSR_plot_solv.jpeg']);

%% Plot r values together
figure(index_n + 2)
tiledlayout(1, 3,'TileSpacing', 'tight','Padding','Tight')
bin_factor = 5;
    
nexttile
plot(N_nums.', r_array(:, 2).', 'Linestyle','-','Marker','.')
hold on
plot(N_nums.', r_array(:, 3).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 4).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 5).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 6).', 'Linestyle','-','Marker','.')
legend('A1g','A2g','B1g','B2g','Eg')
title('Unsolvable wavefunctions')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53', 'HandleVisibility','off')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums.', r_array(:, 7).', 'Linestyle','-','Marker','.')
hold on
plot(N_nums.', r_array(:, 8).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 9).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 10).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 11).', 'Linestyle','-','Marker','.')
legend('A1u','A2u','B1u','B2u','Eu')
title('Unsolvable wavefunctions')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53', 'HandleVisibility','off')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 12).', 'Linestyle','-','Marker','.')
hold on
plot(N_nums.', r_array(:, 13).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 14).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 15).', 'Linestyle','-','Marker','.')
if (l1_even == false && l2_even)
    legend('A1g', 'A1u', 'B2g', 'B2u')
elseif (l1_even && l2_even == false)
    legend('A1g', 'A1u', 'A2g', 'A2u')
elseif (l1_even == false && l2_even == false)
    legend('A1g', 'A1u', 'B1g', 'B1u')
end

title('Solvable wavefunctions')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39', 'HandleVisibility','off')
ylim([ylow_solv, yhigh_solv])
set(figure(index_n + 2),'position',[0,100,1500,300])
saveas(gcf, [folderName '/LSR_plot_combined.jpeg']);

%% Plot proportions of solvable states
figure(index_n + 3)
proportions = sum(solvable_prop, 2);
plot(N_nums(3:end), proportions(3:end), 'Linestyle','-','Marker','.')
title('Proportion of solvable states')
xlabel('N')
prop_solvable = 1 / (4 * (l1 * l1 + 2 * l1 * l2));
prop_range = 0.05;
yline(prop_solvable, '--', ['1/' num2str(4 * (l1 * l1 + 2 * l1 * l2))])
ylim([prop_solvable - prop_range, prop_solvable + prop_range])
ylabel('Proportion')
set(figure(index_n + 3),'position',[0,100,600,400])
saveas(gcf, [folderName '/solvable_proportion.jpeg']);

%% Save variables
save([folderName '/r_array.mat'],'r_array','-v7.3')
save([folderName '/size_array.mat'],'size_array','-v7.3')
save([folderName '/elevels_a1g.mat'],'elevels_a1g','-v7.3')
save([folderName '/elevels_a1g_solv_sorted.mat'],'elevels_a1g_solv_sorted','-v7.3')
save([folderName '/elevels_a1u.mat'],'elevels_a1u','-v7.3')
save([folderName '/elevels_a1u_solv_sorted.mat'],'elevels_a1u_solv_sorted','-v7.3')
save([folderName '/elevels_a2g.mat'],'elevels_a2g','-v7.3')
save([folderName '/elevels_a2u.mat'],'elevels_a2u','-v7.3')
save([folderName '/elevels_b1g.mat'],'elevels_b1g','-v7.3')
save([folderName '/elevels_b1u.mat'],'elevels_b1u','-v7.3')
save([folderName '/elevels_b2g.mat'],'elevels_b2g','-v7.3')
save([folderName '/elevels_b2u.mat'],'elevels_b2u','-v7.3')
save([folderName '/elevels_eg.mat'],'elevels_eg','-v7.3')
save([folderName '/elevels_eu.mat'],'elevels_eu','-v7.3')

if (l1_even == false && l2_even)
    save([folderName '/elevels_b2g_solv_sorted.mat'],'elevels_b2g_solv_sorted','-v7.3')
    save([folderName '/elevels_b2u_solv_sorted.mat'],'elevels_b2u_solv_sorted','-v7.3')
elseif (l1_even && l2_even == false)
    save([folderName '/elevels_a2u_solv_sorted.mat'],'elevels_a2u_solv_sorted','-v7.3')
    save([folderName '/elevels_a2g_solv_sorted.mat'],'elevels_a2g_solv_sorted','-v7.3')
elseif (l1_even == false && l2_even == false)
    save([folderName '/elevels_b1u_solv_sorted.mat'],'elevels_b1u_solv_sorted','-v7.3')
    save([folderName '/elevels_b1g_solv_sorted.mat'],'elevels_b1g_solv_sorted','-v7.3')
end


%% Functions
% Gives x coordinate of index within each face
function x = xfacecoord(index, N, l1, l2, total_sites)
    l3 = l1;
    if (index <= l1 * (l2 + l3) * N^2) % In face 1 or 2
        x = mod(index - 1, l1 * N) + 1;
    elseif (index > total_sites - l1 * l2 * N^2) % In face 6
        x = mod(index - (total_sites - l1 * l2 * N^2) - 1, l1 * N) + 1;
    elseif (mod(index - l1 * (l2 + l3) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % In face 3
        x = mod(index - l1 * (l2 + l3) * N^2 - 1, (2 * l2 + l1) * N) + 1;
    elseif (mod(index - l1 * (l2 + l3) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % In face 4
        x = mod(index - l1 * (l2 + l3) * N^2 - 1, (2 * l2 + l1) * N) + 1 - l2 * N;
    else % In face 5
        x = mod(index - l1 * (l2 + l3) * N^2 - 1, (2 * l2 + l1) * N) + 1 - (l2 + l1) * N;
    end 
end

% Gives y coordinate of index within each face
function y = yfacecoord(index, N, l1, l2, total_sites) 
    l3 = l1;
    if (index <= l1 * l3 * N^2) % in face 1
        y = ceil(index / (l1 * N));
    elseif (index <= l1 * (l2 + l3) * N^2) % in face 2 
        y = ceil((index - l1 * l3 * N^2) / (l1 * N));
    elseif (index <= total_sites - l1 * l2 * N^2) % in faces 3, 4, or 5
        y = ceil((index - l1 * (l2 + l3) * N^2) / ((2 * l2 + l1) * N));
    else % in face 6
        y = ceil((index - l1 * (l2 + l3) * N^2 - l3 * (2 * l2 + l1) * N^2) / (l1 * N));
    end
end

% From face coordinate and face number give index
function i = index_from_coord(x, y, face, N, l1, l2, total_sites)
    l3 = l1;
    if (face == 1)
        i = x + (l1 * N) * (y - 1);
    elseif (face == 2)
        i = (l1 * l3 * N^2) + x + (l1 * N) * (y - 1);
    elseif (face == 3)
        i = l1 * (l3 + l2) * N^2 + x + (2 * l2 + l1) * N * (y - 1);
    elseif (face == 4)
        i = l1 * (l3 + l2) * N^2 + l2 * N + x + (2 * l2 + l1) * N * (y - 1);
    elseif (face == 5)
        i = l1 * (l3 + l2) * N^2 + (l2 + l1) * N + x + (2 * l2 + l1) * N * (y - 1);
    else 
        i = total_sites - l1 * l2 * N^2 + x + l1 * N * (y - 1);
    end
end

% Given index on discretized cube returns new index after sigma_d symmetry
function i = sigma_d_index(index, N, l1, l2, total_sites)
    if (index <= l1 * l1 * N^2) % in face 1
        new_face = 1;
        new_x = l1 * N + 1 - yfacecoord(index, N, l1, l2, total_sites);
        new_y = l1 * N + 1 - xfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (index <= l1 * (l1 + l2) * N^2) % in face 2 
        new_face = 3;
        new_x = yfacecoord(index, N, l1, l2, total_sites);
        new_y = xfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (index > total_sites - l1 * l2 * N^2) % in face 6
        new_face = 5;
        new_x = yfacecoord(index, N, l1, l2, total_sites);
        new_y = xfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (mod(index - l1 * (l1 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % in face 3
        new_face = 2;
        new_x = yfacecoord(index, N, l1, l2, total_sites);
        new_y = xfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (mod(index - l1 * (l1 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % in face 4
        new_face = 4;
        new_x = yfacecoord(index, N, l1, l2, total_sites);
        new_y = xfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    else  % in face 5
        new_face = 6;
        new_x = yfacecoord(index, N, l1, l2, total_sites);
        new_y = xfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    end
end

% Given index on discretized cube returns new index after C2x symmetry
function i = c2_oo_index(index, N, l1, l2, total_sites)
    if (index <= l1 * l1 * N^2) % in face 1
        new_face = 4;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, total_sites);
        new_y = l1 * N + 1 - yfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (index <= l1 * (l1 + l2) * N^2) % in face 2 
        new_face = 2;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, total_sites);
        new_y = l2 * N + 1 - yfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (index > total_sites - l1 * l2 * N^2) % in face 6
        new_face = 6;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, total_sites);
        new_y = l2 * N + 1 - yfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (mod(index - l1 * (l1 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % in face 3
        new_face = 5;
        new_x = xfacecoord(index, N, l1, l2, total_sites);
        new_y = yfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (mod(index - l1 * (l1 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % in face 4
        new_face = 1;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, total_sites);
        new_y = l1 * N + 1 - yfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    else  % in face 5
        new_face = 3;
        new_x = xfacecoord(index, N, l1, l2, total_sites);
        new_y = yfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    end
end

% Given index on discretized cube returns new index after sigma_v symmetry
function i = inv_index(index, N, l1, l2, total_sites)
    if (index <= l1 * l1 * N^2) % in face 1
        new_face = 4;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, total_sites);
        new_y = yfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (index <= l1 * (l1 + l2) * N^2) % in face 2 
        new_face = 6;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, total_sites);
        new_y = yfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (index > total_sites - l1 * l2 * N^2) % in face 6
        new_face = 2;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, total_sites);
        new_y = yfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (mod(index - l1 * (l1 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % in face 3
        new_face = 5;
        new_x = xfacecoord(index, N, l1, l2, total_sites);
        new_y = l1 * N + 1 - yfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (mod(index - l1 * (l1 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % in face 4
        new_face = 1;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, total_sites);
        new_y = yfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    else  % in face 5
        new_face = 3;
        new_x = xfacecoord(index, N, l1, l2, total_sites);
        new_y = l1 * N + 1 - yfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    end
end

% Given index on discretized cube returns new index after C2y symmetry
function i = c2_o_index(index, N, l1, l2, total_sites)
    if (index <= l1 * l1 * N^2) % in face 1
        new_face = 4;
        new_x = yfacecoord(index, N, l1, l2, total_sites);
        new_y = l1 * N + 1 - xfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (index <= l1 * (l1 + l2) * N^2) % in face 2 
        new_face = 5;
        new_x = yfacecoord(index, N, l1, l2, total_sites);
        new_y = l1 * N + 1 - xfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (index > total_sites - l1 * l2 * N^2) % in face 6
        new_face = 3;
        new_x = yfacecoord(index, N, l1, l2, total_sites);
        new_y = l1 * N + 1 - xfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (mod(index - l1 * (l1 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % in face 3
        new_face = 6;
        new_x = l1 * N + 1 - yfacecoord(index, N, l1, l2, total_sites);
        new_y = xfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    elseif (mod(index - l1 * (l1 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % in face 4
        new_face = 1;
        new_x = l1 * N + 1 - yfacecoord(index, N, l1, l2, total_sites);
        new_y = xfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    else  % in face 5
        new_face = 2;
        new_x = l1 * N + 1 - yfacecoord(index, N, l1, l2, total_sites);
        new_y = xfacecoord(index, N, l1, l2, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, total_sites);
    end
end

% Functions to find indices of H matrix
function y = upper(x, N, l1, l2, total_sites) 
    l3 = l1;
    if (x > l1 * (l2 + l3) * N^2 - l1 * N && x <= l1 * (l2 + l3) * N^2) % x in top row of lower arm 
        y = l1 * (l2 + l3) * N^2 + l2 * N + mod(x - 1, l1 * N) + 1;
    elseif (x > l1 * (l2 + l3) * N^2 && x <= (total_sites - l2 * l1 * N^2 - (2 * l2 + l1) * N)) % x in middle arm but not top row
        y = x + (2 * l2 + l1) * N;
    elseif (x > total_sites - l1 * l2 * N^2 - (2 * l2 + l1) * N && x <= total_sites - l1 * l2 * N^2 - (l2 + l1) * N) % x top row of middle arms, left side
        dist_from_vertex = x - (total_sites - l1 * l2 * N^2 - (2 * l2 + l1) * N); % distance from leftmost vertex
        y = total_sites - dist_from_vertex * l1 * N + 1; 
    elseif (x > total_sites - l1 * l2 * N^2 - (2 * l2 + l1) * N && x <= total_sites - l1 * l2 * N^2 - l2 * N) % x in top row of middle arms, center side
        y = x + (l1 + l2) * N;
    elseif (x > total_sites - l1 * l2 * N^2 - (2 * l2 + l1) * N && x <= total_sites - l1 * l2 * N^2) % x in top row of middle arms, right side
        dist_from_fold = x - (total_sites - l1 * l2 * N^2 - l2 * N);
        y = total_sites - l1 * l2 * N^2 + dist_from_fold * l1 * N;
    elseif (x > total_sites - l1 * N) % x in top row of upper arm
        y = mod((x - total_sites - l1 * N) - 1, l1 * N) + 1;
    else % x in non-top row upper arm or non-top row lower arm
        y = x + l1 * N;
    end
end

function y = lower(x, N, l1, l2, total_sites) 
    l3 = l1;
    if (x > 0 && x <= l1 * N) % x in lower arm bottom row
        y = total_sites - l1 * N + mod(x - 1, l1 * N) + 1;
    elseif (x > l1 * (l3 + l2) * N^2 && x <= l1 * (l3 + l2) * N^2 + l2 * N) % x in middle arms, bottom row, left side
        dist_from_vertex = x - l1 * (l3 + l2) * N^2;
        y = l1 * l3 * N^2 + (dist_from_vertex - 1) * l1 * N + 1;
    elseif (x > l1 * (l3 + l2) * N^2 && x <= l1 * (l3 + l2) * N^2 + (l2 + l1) * N) % x in middle arms, bottom row, center
        y = x - (l2 + l1) * N;
    elseif (x > l1 * (l3 + l2) * N^2 && x <= l1 * (l3 + l2) * N^2 + (2 * l2 + l1) * N) % x in middle arms, bottom row, right side
        dist_from_fold = x - (l1 * (l3 + l2) * N^2 + (l2 + l1) * N);
        y = l1 * (l3 + l2) * N^2 - (dist_from_fold - 1) * l1 * N ;
    elseif (x > l1 * (l3 + l2) * N^2 + (2 * l2 + l1) * N && x <= total_sites - l1 * l2 * N^2) % x in middle arms, above bottom row
        y = x - (2 * l2 + l1) * N;
    elseif (x > total_sites - l1 * l2 * N^2 && x <= total_sites - l1 * l2 * N^2 + l1 * N) % x in upper arm, bottom row
        y = x - (l1 + l2) * N;
    else % x in upper arm, above bottom row or in lower arm above bottom row
        y = x - l1 * N;
    end
end

function y = right(x, N, l1, l2, total_sites)
    l3 = l1;
    if (x <= l3 * l1 * N^2 && floor(x / (l1 * N)) == x / (l1 * N)) % x in lower half of lower arm, rightmost column
        row_from_bottom = x / (l1 * N);
        y = total_sites - l2 * l1 * N^2 - (row_from_bottom - 1) * (2 * l2 + l1) * N;
    elseif (x > l1 * l3 * N^2 && x <= l1 * (l2 + l3) * N^2 && floor(x / (l1 * N)) == x / (l1 * N)) % x in upper half of lower arm, rightmost column
        row_from_top = (l1 * (l2 + l3) * N^2 - x) / (l1 * N) + 1;
        y = (l3 + l2) * l1 * N^2 + (l2 + l1) * N + row_from_top;
    elseif (x > (l3 + l2) * l1 * N^2 && x <= total_sites - l2 * l1 * N^2 && floor((x - (l2 + l3) * l1 * N^2)/((2 * l2 + l1) * N)) == (x - (l2 + l3) * l1 * N^2)/((2 * l2 + l1) * N)) % x in middle arm, rightmost column 
        row_from_top_2 = (total_sites - l2 * l1 * N^2 - x)/ ((2 * l2 + l1) * N) + 1;
        y = row_from_top_2 * l1 * N;
    elseif (x > total_sites - l2 * l1 * N^2 && floor((x - (total_sites - l2 * l1 * N^2)) / (l1 * N)) == (x - (total_sites - l2 * l1 * N^2)) / (l1 * N))  % x in upper arm, rightmost column 
        row_from_bottom_2 = (x - (total_sites - l2 * l1 * N^2)) / (l1 * N);
        y = total_sites - l2 * l1 * N^2 - l2 * N + row_from_bottom_2;
    else 
        y = x + 1;
    end
end

function y = left(x, N, l1, l2, total_sites)
    l3 = l1;
    if (x <= l3 * l1 * N^2 && (floor((x - 1) / (l1 * N)) == (x - 1) / (l1 * N))) % x in lower half of lower arm, leftmost column
        row_from_bottom = (x - 1) / (l1 * N) + 1;
        y = total_sites - l2 * l1 * N^2 - (row_from_bottom) * (2 * l2 + l1) * N + 1;
    elseif (x > l3 * l1 * N^2 && x <= (l3 + l2) * l1 * N^2 && floor((x - 1) / (l1 * N)) == (x - 1) / (l1 * N)) % x in upper half of lower arm, leftmost column
        row_from_bottom_2 = ceil(x / (l1 * N)) - l3 * N;
        y = (l3 + l2) * l1 * N^2 + row_from_bottom_2;
    elseif (x > (l3 + l2) * l1 * N^2 && x <= total_sites - l2 * l1 * N^2 && floor((x - 1 - (l3 + l2) * l1 * N^2)/((2 * l2 + l1) * N)) == (x - 1 - (l2 + l3) * l1 * N^2)/((2 * l2 + l1) * N)) % x in middle arm, leftmost column 
        row_from_top = ceil((total_sites - l2 * l1 * N^2 - x)/((2 * l2 + l1) * N));
        y = l1 * N * (row_from_top - 1) + 1;
    elseif (x > total_sites - l2 * l1 * N^2 && floor((x - (total_sites - l2 * l1 * N^2) - 1) / (l1 * N)) == (x - (total_sites - l2 * l1 * N^2) - 1) / (l1 * N))  % x in upper arm, leftmost column 
        row_from_top_2 = ceil((total_sites - x) / (l1 * N));
        y = total_sites - l2 * l1 * N^2 - (2 * l2 + l1) * N + row_from_top_2;
    else  
        y = x - 1;
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