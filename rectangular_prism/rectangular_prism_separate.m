%% Calculate r values for various N
clear
close all
tolerance = 1e-8;  % Tolerance for degenerate energy levels
corner_potential = -0;   % Corner potential

% Side length ratios
l1 = 1;
l2 = 2;
l3 = 3;

% length parity indicators
all_same_parity = false;
l1_l2_same_parity = false;
l1_l3_same_parity = false;
l2_l3_same_parity = false;

l1_l2_even = false;
l1_l3_even = false;
l2_l3_even = false;

l1_even = false;
l2_even = false;
l3_even = false;

if (mod(l1, 2) == 0)
    l1_even = true;
elseif (mod(l2, 2) == 0)
    l2_even = true;
elseif (mod(l3, 2) == 0)
    l3_even = true;
end

if (mod(l1 - l2, 2) == 0) && (mod(l1 - l3, 2) == 0)
    all_same_parity = true;
elseif (mod(l1 - l2, 2) == 0)
    l1_l2_same_parity = true;
    if (l1_even)
        l1_l2_even = true;
    end
elseif (mod(l1 - l3, 2) == 0)
    l1_l3_same_parity = true;
    if (mod(l1, 2) == 0)
        l1_l3_even = true;
    end
else
    l2_l3_same_parity = true;
    if (l2_even)
        l2_l3_even = true;
    end
end


% Find max N
upper_lim = 30000; % Max number of sites allowed in memory
max_N = floor(sqrt(upper_lim / (2 * ((l1 * l2) + (l2 * l3) + (l3 * l1)))));
fprintf(['Max N = ' num2str(max_N) '\n'])
N_nums = primes(max_N);
%N_nums = [2];

% Remove small primes (2, 3, 5)
N_nums = N_nums(4:end);

% Make directory
folderName = ['l1=' num2str(l1) '_l2=' num2str(l2) '_l3=' num2str(l3) '_separated'];
mkdir(folderName);

% Allocate variable space
if (all_same_parity == true)
    r_array = zeros(size(N_nums, 1), 11);
    size_array = zeros(size(N_nums, 1), 11);
    solvable_prop = zeros(size(N_nums, 1), 2);
else 
    r_array = zeros(size(N_nums, 1), 13);
    size_array = zeros(size(N_nums, 1), 13);
    solvable_prop = zeros(size(N_nums, 1), 4);
end

%% Begin Program
% Set indices
index_n = 1;
index_r = 1;
index_size = 1;

%% Start loop
for n = 1:size(N_nums, 2)

    N = N_nums(n);
    fprintf(['\n    N = ' num2str(N) '\n'])
    total_sites = 2 * ((l1 * l2) + (l2 * l3) + (l3 * l1)) * N^2; % Total number of sites in unfolded polyhedron

    %% Diagonalize H and sigma_d matrices
    tic;

    % Generate Hamiltonian matrix
    L = 1/N; % Spacing
    H = zeros(total_sites, total_sites); % Hamiltonian matrix

    %for i = 1: total_sites
        %fprintf([num2str(i) ' ' num2str(upper(i, N, l1, l2, l3, total_sites)) ' ' num2str(lower(i, N, l1, l2, l3, total_sites)) ' ' ...
            %num2str(right(i, N, l1, l2, l3, total_sites)) ' ' num2str(left(i, N, l1, l2, l3, total_sites)) '\n'])
    %end
    
    for i = 1:total_sites
        H(i, i) = 4/L^2;
        H(i, upper(i, N, l1, l2, l3, total_sites)) = -1/L^2;
        H(i, lower(i, N, l1, l2, l3, total_sites)) = -1/L^2;
        H(i, right(i, N, l1, l2, l3, total_sites)) = -1/L^2;
        H(i, left(i, N, l1, l2, l3, total_sites)) = -1/L^2;       
    end
    
    %% Generate matrix for c2x symmetry 
    C2x = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        C2x(i, c2x_index(i, N, l1, l2, l3, total_sites)) = 1;
        %fprintf([num2str(i) ' ' num2str(c2x_index(i, N, l1, l2, l3, total_sites)) '\n'])
    end
    
    %% Calculate eigenvectors and eigenvalues
    [eigenvectors, eigenvalues] = eig(H + ((1 + sqrt(5))/2) * C2x);
    fprintf('evals and evecs done \n')
    toc
    
    %% Generate Symmetry Matrices
    
    % Generate matrix for C2y
    C2y = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        C2y(i, c2y_index(i, N, l1, l2, l3, total_sites)) = 1;
        %fprintf([num2str(i) ' ' num2str(c2y_index(i, N, l1, l2, l3, total_sites)) '\n'])
    end

    % Generate matrix for C2z
    C2z = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        C2z(i, c2z_index(i, N, l1, l2, l3, total_sites)) = 1;
        %fprintf([num2str(i) ' ' num2str(c2z_index(i, N, l1, l2, l3, total_sites)) '\n'])
    end

    % Generate matrix for i
    inv = zeros(total_sites, total_sites);
    for i = 1:(total_sites)
        inv(i, inv_index(i, N, l1, l2, l3, total_sites)) = 1;
        %fprintf([num2str(i) ' ' num2str(inv_index(i, N, l1, l2, l3, total_sites)) '\n'])
    end

    
    %% Calculate characters of symmetries
    tic;
    % Preallocate data array
    c2x_evals = zeros(total_sites, 1);
    c2y_evals = zeros(total_sites, 1);
    c2z_evals = zeros(total_sites, 1);
    inv_evals = zeros(total_sites, 1);
    
    % Find energy and c2x eigenvluaes
    for i = 1:(total_sites)
        c2x_evals(i) = (eigenvectors(:, i).' * C2x * eigenvectors(:, i));
        eigenvalues(i, i) = eigenvalues(i, i) - ((1 + sqrt(5))/2) * c2x_evals(i);
    end
    fprintf('c2x done \n')
    
    % Find c2y eigenvalues
    for i = 1:(total_sites)
        c2y_evals(i) = (eigenvectors(:, i).' * C2y * eigenvectors(:, i));
    end
    fprintf('c2y done \n')

    % Find c2z eigenvalues
    for i = 1:(total_sites)
        c2z_evals(i) = (eigenvectors(:, i).' * C2z * eigenvectors(:, i));
    end
    fprintf('c2z done \n')

    % Find inv eigenvalues
    for i = 1:(total_sites)
        inv_evals(i) = (eigenvectors(:, i).' * inv * eigenvectors(:, i));
    end
    fprintf('inv done \n')
    
    
    %% Analyze Degeneracies
    % Reorder energies and evecs
    energy_levels = zeros(total_sites, 7);
    [energies, id] = sort(diag(eigenvalues));
    eigenvectors = eigenvectors(:, id);
    c2x_evals_sorted = c2x_evals(id);
    c2y_evals_sorted = c2y_evals(id);
    c2z_evals_sorted = c2z_evals(id);
    inv_evals_sorted = inv_evals(id);

    % Analyze for degeneracies
    energy = energies(1);
    degeneracy = 1;
    index = 1;
    trace_c2x = c2x_evals_sorted(1);
    trace_c2y = c2y_evals_sorted(1);
    trace_c2z = c2z_evals_sorted(1);
    trace_inv = inv_evals_sorted(1);

    for i = 2:(total_sites)
        if (abs(energies(i) - energy) < tolerance) % Next energy is degenerate
            degeneracy = degeneracy + 1;
            trace_c2x = trace_c2x + c2x_evals_sorted(i);
            trace_c2y = trace_c2y + c2y_evals_sorted(i);
            trace_c2z = trace_c2z + c2z_evals_sorted(i);
            trace_inv = trace_inv + inv_evals_sorted(i);
        else % Next energy is new
            % Record stats for previous energy level
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_c2x;
            energy_levels(index, 4) = trace_c2y;
            energy_levels(index, 5) = trace_c2z;
            energy_levels(index, 6) = trace_inv;
            energy_levels(index, 7) = i - 1;
    
            % Reset variables for new energy level
            energy = energies(i);
            degeneracy = 1;
            index = index + 1;
            trace_c2x = c2x_evals_sorted(i);
            trace_c2y = c2y_evals_sorted(i);
            trace_c2z = c2z_evals_sorted(i);
            trace_inv = inv_evals_sorted(i);            
        end
    
        % Record energy level if we reach the end
        if (i == total_sites)
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_c2x;
            energy_levels(index, 4) = trace_c2y;
            energy_levels(index, 5) = trace_c2z;
            energy_levels(index, 6) = trace_inv;
            energy_levels(index, 7) = i;
            
        end
    end
    
    energy_levels = energy_levels(1:index, :);
    
    %% Separate the representation classes
    % Allocate space for each irrep
    elevels_ag = zeros(total_sites, 7);
    elevels_b1g = zeros(total_sites, 7);
    elevels_b2g = zeros(total_sites, 7);
    elevels_b3g = zeros(total_sites, 7);
    elevels_au = zeros(total_sites, 7);
    elevels_b1u = zeros(total_sites, 7);
    elevels_b2u = zeros(total_sites, 7);
    elevels_b3u = zeros(total_sites, 7);

    % Allocate space for solvable elevels
    elevels_ag_solv = zeros(total_sites, 7);
    elevels_au_solv = zeros(total_sites, 7);

    if (l1_l2_same_parity == true)
        elevels_b3g_solv = zeros(total_sites, 7);
        elevels_b3u_solv = zeros(total_sites, 7);
    elseif (l1_l3_same_parity == true)
        elevels_b2g_solv = zeros(total_sites, 7);
        elevels_b2u_solv = zeros(total_sites, 7);
    elseif (l2_l3_same_parity == true)
        elevels_b1g_solv = zeros(total_sites, 7);
        elevels_b1u_solv = zeros(total_sites, 7);
    end

    % Set indices
    index_ag = 1;
    index_b1g = 1;
    index_b2g = 1;
    index_b3g = 1;
    index_au = 1;
    index_b1u = 1;
    index_b2u = 1;
    index_b3u = 1;

    index_ag_solv = 1;
    index_au_solv = 1;

    if (l1_l2_same_parity == true)
        index_b3g_solv = 1;
        index_b3u_solv = 1;
    elseif (l1_l3_same_parity == true)
        index_b2g_solv = 1;
        index_b2u_solv = 1;
    elseif (l2_l3_same_parity == true)
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
        trace_c2x = energy_levels_rounded(i, 3);
        trace_c2y = energy_levels_rounded(i, 4);
        trace_c2z = energy_levels_rounded(i, 5);
        trace_inv = energy_levels_rounded(i, 6);
    
        traces = [trace_e, trace_c2x, trace_c2y, trace_c2z, trace_inv];
    
        if (isequal(traces, [1, 1, 1, 1, 1])) % Ag
            elevels_ag(index_ag, :) = energy_levels_rounded(i, :);
            index_ag = index_ag + 1;
        elseif (isequal(traces, [1, -1, -1, 1, 1])) % B1g
            elevels_b1g(index_b1g, :) = energy_levels_rounded(i, :);
            index_b1g = index_b1g + 1;
        elseif (isequal(traces, [1, -1, 1, -1, 1])) % B2g
            elevels_b2g(index_b2g, :) = energy_levels_rounded(i, :);
            index_b2g = index_b2g + 1;
        elseif (isequal(traces, [1, 1, -1, -1, 1])) % B3g
            elevels_b3g(index_b3g, :) = energy_levels_rounded(i, :);
            index_b3g = index_b3g + 1;
        elseif (isequal(traces, [1, 1, 1, 1, -1])) % Au
            elevels_au(index_au, :) = energy_levels_rounded(i, :);
            index_au = index_au + 1;
        elseif (isequal(traces, [1, -1, -1, 1, -1])) % B1u
            elevels_b1u(index_b1u, :) = energy_levels_rounded(i, :);
            index_b1u = index_b1u + 1;
        elseif (isequal(traces, [1, -1, 1, -1, -1])) % B2u
            elevels_b2u(index_b2u, :) = energy_levels_rounded(i, :);
            index_b2u = index_b2u + 1;
        elseif (isequal(traces, [1, 1, -1, -1, -1])) % B3u
            elevels_b3u(index_b3u, :) = energy_levels_rounded(i, :);
            index_b3u = index_b3u + 1;
        elseif (isequal(traces, [2, 2, 2, 2, 0])) % Accidental degeneracy Ag + Au
            elevels_ag_solv(index_ag_solv, :) = energy_levels_rounded(i, :);
            index_ag_solv = index_ag_solv + 1;
            elevels_au_solv(index_au_solv, :) = energy_levels_rounded(i, :);
            index_au_solv = index_au_solv + 1;
        elseif (isequal(traces, [2, -2, -2, 2, 0])) % Accidental degeneracy B1g + B1u
            elevels_b1g_solv(index_b1g_solv, :) = energy_levels_rounded(i, :);
            index_b1g_solv = index_b1g_solv + 1;
            elevels_b1u_solv(index_b1u_solv, :) = energy_levels_rounded(i, :);
            index_b1u_solv = index_b1u_solv + 1;
        elseif (isequal(traces, [2, -2, 2, -2, 0])) % Accidental degeneracy B2g + B2u
            elevels_b2g_solv(index_b2g_solv, :) = energy_levels_rounded(i, :);
            index_b2g_solv = index_b2g_solv + 1;
            elevels_b2u_solv(index_b2u_solv, :) = energy_levels_rounded(i, :);
            index_b2u_solv = index_b2u_solv + 1;
        elseif (isequal(traces, [2, 2, -2, -2, 0])) % Accidental degeneracy B3g + B3u
            elevels_b3g_solv(index_b3g_solv, :) = energy_levels_rounded(i, :);
            index_b3g_solv = index_b3g_solv + 1;
            elevels_b3u_solv(index_b3u_solv, :) = energy_levels_rounded(i, :);
            index_b3u_solv = index_b3u_solv + 1;
        else 
            fprintf([num2str(energy_levels_rounded(i, 1)) ' ' num2str(energy_levels_rounded(i, 2)) ' ' ...
                num2str(energy_levels_rounded(i, 3)) ' ' num2str(energy_levels_rounded(i, 4)) ' ' ...
                num2str(energy_levels_rounded(i, 5)) ' ' num2str(energy_levels_rounded(i, 6)) ' ' num2str(energy_levels_rounded(i, 7)) '\n'])
        end
    end

    %% Calculate energies of non-degen solvable wavefunctions
    % Calculate solvable wavefunctions energies from Ag (n1, n1) with n1 even
    ag_exp_energies = zeros(ceil(N / 2), 1);
    for i = 1:ceil(N / 2)
        n1 = (i - 1) * 2;
        ag_exp_energies(i) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
    end

    % Calculate solvable wavefunction energies (n1, n1) with n1 odd
    if (all_same_parity == true) % Energies in Au (all same parity)
        au_exp_energies = zeros(ceil(N / 2), 1);
        for i = 1:ceil(N / 2)
            n1 = (i - 1) * 2;
            au_exp_energies(i) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
        end
    elseif (l1_l2_same_parity == true) % Energies in B3 (l1 and l2 same parity)
        if (l1_l2_even == false) % B3g
            b3g_exp_energies = zeros(ceil(N / 2), 1);
            for i = 1:ceil(N / 2)
                n1 = (i - 1) * 2 + 1;
                b3g_exp_energies(i) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
            end
        else % B3u
            b3u_exp_energies = zeros(ceil(N / 2), 1);
            for i = 1:ceil(N / 2)
                n1 = (i - 1) * 2 + 1;
                b3u_exp_energies(i) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
            end
        end
    elseif (l1_l3_same_parity == true) % Energies in B2 (l1 and l3 same parity)
        if (l1_l3_even == false) % In B2g
            b2g_exp_energies = zeros(ceil(N / 2), 1);
            for i = 1:ceil(N / 2)
                n1 = (i - 1) * 2 + 1;
                b2g_exp_energies(i) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
            end
        else % In B2u
            b2u_exp_energies = zeros(ceil(N / 2), 1);
            for i = 1:ceil(N / 2)
                n1 = (i - 1) * 2 + 1;
                b2u_exp_energies(i) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
            end
        end
    elseif (l2_l3_same_parity == false) % Energies in B1 (l2 and l3 same parity)
        if (l2_l3_even == true) % In B1g
            b1g_exp_energies = zeros(ceil(N / 2), 1);
            for i = 1:ceil(N / 2)
                n1 = (i - 1) * 2 + 1;
                b1g_exp_energies(i) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
            end
        else % In B1u
            b1u_exp_energies = zeros(ceil(N / 2), 1);
            for i = 1:ceil(N / 2)
                n1 = (i - 1) * 2 + 1;
                b1u_exp_energies(i) = 2 * N^2 * (2 - cos(pi * n1 / N) - cos(pi * n1 / N));
            end
        end
    end
    
    %% Remove non-degen solvable wavefunctions

    % Remove from Ag
    for i = 1:(index_ag - 1)
        if (i < index_ag)
            for j = 1:ceil(N / 2)
                if (abs(ag_exp_energies(j) - elevels_ag(i, 1)) < tolerance)
                    % Add to elevel_ag_solv
                    elevels_ag_solv(index_ag_solv, :) = elevels_ag(i, :);
                    index_ag_solv = index_ag_solv + 1;
    
                    % Delete from elevels_ag
                    elevels_ag(i, :) = [];
                    index_ag = index_ag - 1;
                end
            end 
        end
    end

    if (all_same_parity == true) % Remove from Au (even n1)
        for i = 1:(index_au - 1)
            if (i < index_au)
                for j = 1:ceil(N / 2)
                    if (abs(au_exp_energies(j) - elevels_au(i, 1)) < tolerance)
                        % Add to elevel_au_solv
                        elevels_au_solv(index_au_solv, :) = elevels_au(i, :);
                        index_au_solv = index_au_solv + 1;
    
                        % Delete from elevels_au
                        elevels_au(i, :) = [];
                        index_au = index_au - 1;
                    end
                end
            end
        end
    elseif (l1_l2_same_parity == true) % Remove from B irreps (odd n1)
        if (l1_l2_even == false) % Remove from B3g
            for i = 1:(index_b3g - 1)
                if (i < index_b3g)
                    for j = 1:ceil(N / 2)
                        if (abs(b3g_exp_energies(j) - elevels_b3g(i, 1)) < tolerance)
                            % Add to elevel_b3g_solv
                            elevels_b3g_solv(index_b3g_solv, :) = elevels_b3g(i, :);
                            index_b3g_solv = index_b3g_solv + 1;
        
                            % Delete from elevels_b3g
                            elevels_b3g(i, :) = [];
                            index_b3g = index_b3g - 1;
                        end
                    end
                end
            end
        else
            for i = 1:(index_b3u - 1) % Remove from B3u
                if (i < index_b3u)
                    for j = 1:ceil(N / 2)
                        if (abs(b3u_exp_energies(j) - elevels_b3u(i, 1)) < tolerance)
                            % Add to elevel_b3u_solv
                            elevels_b3u_solv(index_b3u_solv, :) = elevels_b3u(i, :);
                            index_b3u_solv = index_b3u_solv + 1;
        
                            % Delete from elevels_b3u
                            elevels_b3u(i, :) = [];
                            index_b3u = index_b3u - 1;
                        end
                    end
                end
            end
        end
    elseif (l1_l3_same_parity == true) % Remove from B2 (odd n1)
        if (l1_l3_even == false) % Remove from B2g
            for i = 1:(index_b2g - 1)
                if (i < index_b2g)
                    for j = 1:ceil(N / 2)
                        if (abs(b2g_exp_energies(j) - elevels_b2g(i, 1)) < tolerance)
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
        else
            for i = 1:(index_b2u - 1) % Remove from B2u
                if (i < index_b2u)
                    for j = 1:ceil(N / 2)
                        if (abs(b2u_exp_energies(j) - elevels_b2u(i, 1)) < tolerance)
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
        end
    elseif (l2_l3_same_parity == true) % Remove from B1
        if (l2_l3_even == false) % Remove from B1g
            for i = 1:(index_b1g - 1)
                if (i < index_b1g)
                    for j = 1:ceil(N / 2)
                        if (abs(b1g_exp_energies(j) - elevels_b1g(i, 1)) < tolerance)
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
        else
            for i = 1:(index_b1u - 1) % Remove from B1u
                if (i < index_b1u)
                    for j = 1:ceil(N / 2)
                        if (abs(b1u_exp_energies(j) - elevels_b1u(i, 1)) < tolerance)
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
        end
    end

    
    %% Remove extra rows of zeros
    if (index_ag > 1)
        elevels_ag = elevels_ag(1:(index_ag - 1), :);
    end 
    
    if (index_b1g > 1)
        elevels_b1g = elevels_b1g(1:(index_b1g - 1), :);
    end
    
    if (index_b2g > 1)
        elevels_b2g = elevels_b2g(1:(index_b2g - 1), :);
    end
    
    if (index_b3g > 1)
        elevels_b3g = elevels_b3g(1:(index_b3g - 1), :);
    end
    
    if (index_au > 1)
        elevels_au = elevels_au(1:(index_au - 1), :);
    end 
        
    if (index_b1u > 1)
        elevels_b1u = elevels_b1u(1:(index_b1u - 1), :);
    end
    
    if (index_b2u > 1)
        elevels_b2u = elevels_b2u(1:(index_b2u - 1), :);
    end
    
    if (index_b3u > 1)
        elevels_b3u = elevels_b3u(1:(index_b3u - 1), :);
    end

    if (index_ag_solv > 1)
        elevels_ag_solv = elevels_ag_solv(1:(index_ag_solv - 1), :);
    end 

    if (index_au_solv > 1)
        elevels_au_solv = elevels_au_solv(1:(index_au_solv - 1), :);
    end 
    
    if (l1_l2_same_parity == true)
        if (index_b3g_solv > 1)
            elevels_b3g_solv = elevels_b3g_solv(1:(index_b3g_solv - 1), :);
        end

        if (index_b3u_solv > 1)
            elevels_b3u_solv = elevels_b3u_solv(1:(index_b3u_solv - 1), :);
        end
    elseif (l1_l3_same_parity == true)
        if (index_b2g_solv > 1)
            elevels_b2g_solv = elevels_b2g_solv(1:(index_b2g_solv - 1), :);
        end

        if (index_b2u_solv > 1)
            elevels_b2u_solv = elevels_b2u_solv(1:(index_b2u_solv - 1), :);
        end
    elseif (l2_l3_same_parity == true)
        if (index_b1g_solv > 1)
            elevels_b1g_solv = elevels_b1g_solv(1:(index_b1g_solv - 1), :);
        end

        if (index_b1u_solv > 1)
            elevels_b1u_solv = elevels_b1u_solv(1:(index_b1u_solv - 1), :);
        end
    end 

    % Resort elevels_solv for Ag, B irreps
    elevels_ag_solv_sorted = sortrows(elevels_ag_solv);

    if (all_same_parity == true) 
        elevels_au_solv_sorted = sortrows(elevels_au_solv);
    elseif (l1_l2_same_parity == true)
        elevels_b3g_solv_sorted = sortrows(elevels_b3g_solv);
        elevels_b3u_solv_sorted = sortrows(elevels_b3u_solv);
    elseif (l1_l3_same_parity == true)
        elevels_b2g_solv_sorted = sortrows(elevels_b2g_solv);
        elevels_b2u_solv_sorted = sortrows(elevels_b2u_solv);
    elseif (l2_l3_same_parity == true)
        elevels_b1g_solv_sorted = sortrows(elevels_b1g_solv);
        elevels_b1u_solv_sorted = sortrows(elevels_b1u_solv);
    end
    
    %% Find energy level spacings
    % Allocate variables
    spacings_ag = zeros(size(elevels_ag, 1) - 1, 1);
    spacings_b1g = zeros(size(elevels_b1g, 1) - 1, 1);
    spacings_b2g = zeros(size(elevels_b2g, 1) - 1, 1);
    spacings_b3g = zeros(size(elevels_b3g, 1) - 1, 1);
    spacings_au = zeros(size(elevels_au, 1) - 1, 1);
    spacings_b1u = zeros(size(elevels_b1u, 1) - 1, 1);
    spacings_b2u = zeros(size(elevels_b2u, 1) - 1, 1);
    spacings_b3u = zeros(size(elevels_b3u, 1) - 1, 1);

    spacings_ag_solv = zeros(size(elevels_ag_solv, 1) - 1, 1);
    spacings_au_solv = zeros(size(elevels_au_solv, 1) - 1, 1);

    if (all_same_parity == true) 
        spacings_au_solv = zeros(size(elevels_au_solv, 1) - 1, 1);
    elseif (l1_l2_same_parity == true)
        spacings_b3g_solv = zeros(size(elevels_b3g_solv, 1) - 1, 1);
        spacings_b3u_solv = zeros(size(elevels_b3u_solv, 1) - 1, 1);
    elseif (l1_l3_same_parity == true)
        spacings_b2g_solv = zeros(size(elevels_b2g_solv, 1) - 1, 1);
        spacings_b2u_solv = zeros(size(elevels_b2u_solv, 1) - 1, 1);
    elseif (l2_l3_same_parity == true)
        spacings_b1g_solv = zeros(size(elevels_b1g_solv, 1) - 1, 1);
        spacings_b1u_solv = zeros(size(elevels_b1u_solv, 1) - 1, 1);
    end

    % Compute spacings
    for i = 1:(size(elevels_ag, 1) - 1)
        spacings_ag(i) = abs(elevels_ag(i, 1) - elevels_ag(i+1, 1));
    end
    
    for i = 1:(size(elevels_b1g, 1) - 1)
        spacings_b1g(i) = abs(elevels_b1g(i, 1) - elevels_b1g(i+1, 1));
    end
    
    for i = 1:(size(elevels_b2g, 1) - 1)
        spacings_b2g(i) = abs(elevels_b2g(i, 1) - elevels_b2g(i+1, 1));
    end
    
    for i = 1:(size(elevels_b3g, 1) - 1)
        spacings_b3g(i) = abs(elevels_b3g(i, 1) - elevels_b3g(i+1, 1));
    end
    
    for i = 1:(size(elevels_au, 1) - 1)
        spacings_au(i) = abs(elevels_au(i, 1) - elevels_au(i+1, 1));
    end
    
    for i = 1:(size(elevels_b1u, 1) - 1)
        spacings_b1u(i) = abs(elevels_b1u(i, 1) - elevels_b1u(i+1, 1));
    end
    
    for i = 1:(size(elevels_b2u, 1) - 1)
        spacings_b2u(i) = abs(elevels_b2u(i, 1) - elevels_b2u(i+1, 1));
    end
    
    for i = 1:(size(elevels_b3u, 1) - 1)
        spacings_b3u(i) = abs(elevels_b3u(i, 1) - elevels_b3u(i+1, 1));
    end

    for i = 1:(size(elevels_ag_solv_sorted, 1) - 1)
        spacings_ag_solv(i) = abs(elevels_ag_solv_sorted(i, 1) - elevels_ag_solv_sorted(i + 1, 1));
    end

    for i = 1:(size(elevels_au_solv, 1) - 1)
        spacings_au_solv(i) = abs(elevels_au_solv(i, 1) - elevels_au_solv(i+1, 1));
    end

    if (all_same_parity == true) 
        for i = 1:(size(elevels_au_solv_sorted, 1) - 1)
            spacings_au_solv(i) = abs(elevels_au_solv_sorted(i, 1) - elevels_au_solv_sorted(i+1, 1));
        end
    elseif (l1_l2_same_parity == true)
        for i = 1:(size(elevels_b3g_solv_sorted, 1) - 1)
            spacings_b3g_solv(i) = abs(elevels_b3g_solv_sorted(i, 1) - elevels_b3g_solv_sorted(i+1, 1));
        end

        for i = 1:(size(elevels_b3u_solv_sorted, 1) - 1)
            spacings_b3u_solv(i) = abs(elevels_b3u_solv_sorted(i, 1) - elevels_b3u_solv_sorted(i+1, 1));
        end
    elseif (l1_l3_same_parity == true)
        for i = 1:(size(elevels_b2g_solv_sorted, 1) - 1)
            spacings_b2g_solv(i) = abs(elevels_b2g_solv_sorted(i, 1) - elevels_b2g_solv_sorted(i+1, 1));
        end

        for i = 1:(size(elevels_b2u_solv_sorted, 1) - 1)
            spacings_b2u_solv(i) = abs(elevels_b2u_solv_sorted(i, 1) - elevels_b2u_solv_sorted(i+1, 1));
        end
    elseif (l2_l3_same_parity == true)
        for i = 1:(size(elevels_b1g_solv_sorted, 1) - 1)
            spacings_b1g_solv(i) = abs(elevels_b1g_solv_sorted(i, 1) - elevels_b1g_solv_sorted(i+1, 1));
        end

        for i = 1:(size(elevels_b1u_solv_sorted, 1) - 1)
            spacings_b1u_solv(i) = abs(elevels_b1u_solv_sorted(i, 1) - elevels_b1u_solv_sorted(i+1, 1));
        end
    end

    % Compute level spacing ratios
    r_ag = LSR(spacings_ag);
    r_b1g = LSR(spacings_b1g);
    r_b2g = LSR(spacings_b2g);
    r_b3g = LSR(spacings_b3g);
    r_au = LSR(spacings_au);
    r_b1u = LSR(spacings_b1u);
    r_b2u = LSR(spacings_b2u);
    r_b3u = LSR(spacings_b3u);

    r_ag_solv = LSR(spacings_ag_solv);
    r_au_solv = LSR(spacings_au_solv);

    if (l1_l2_same_parity == true)
        r_b3g_solv = LSR(spacings_b3g_solv);
        r_b3u_solv = LSR(spacings_b3u_solv);
    elseif (l1_l3_same_parity == true)
        r_b2g_solv = LSR(spacings_b2g_solv);
        r_b2u_solv = LSR(spacings_b2u_solv);
    elseif (l2_l3_same_parity == true)
        r_b1g_solv = LSR(spacings_b1g_solv);
        r_b1u_solv = LSR(spacings_b1u_solv);
    end

    if (all_same_parity == true) 
        r_array(index_r, :) = [N, r_ag, r_b1g, r_b2g, r_b3g, r_au, r_b1u, r_b2u, r_b3u, r_ag_solv, r_au_solv];
    elseif (l1_l2_same_parity == true)
        r_array(index_r, :) = [N, r_ag, r_b1g, r_b2g, r_b3g, r_au, r_b1u, r_b2u, r_b3u, r_ag_solv, r_au_solv, r_b3g_solv, r_b3u_solv];
    elseif (l1_l3_same_parity == true)
        r_array(index_r, :) = [N, r_ag, r_b1g, r_b2g, r_b3g, r_au, r_b1u, r_b2u, r_b3u, r_ag_solv, r_au_solv, r_b2g_solv, r_b2u_solv];
    elseif (l2_l3_same_parity == true)
        r_array(index_r, :) = [N, r_ag, r_b1g, r_b2g, r_b3g, r_au, r_b1u, r_b2u, r_b3u, r_ag_solv, r_au_solv, r_b1g_solv, r_b1u_solv];
    end
   
    % Fill in size of irreps
    

    if (all_same_parity == true) 
        size_array(index_size, :) = [N, size(elevels_ag, 1), size(elevels_b1g, 1), ...
        size(elevels_b2g, 1), size(elevels_b3g, 1), size(elevels_au, 1), ...
        size(elevels_b1u, 1), size(elevels_b2u, 1), size(elevels_b3u, 1), ...
        size(elevels_ag_solv, 1), size(elevels_au_solv, 1)];
    elseif (l1_l2_same_parity == true)
        size_array(index_size, :) = [N, size(elevels_ag, 1), size(elevels_b1g, 1), ...
        size(elevels_b2g, 1), size(elevels_b3g, 1), size(elevels_au, 1), ...
        size(elevels_b1u, 1), size(elevels_b2u, 1), size(elevels_b3u, 1), ...
        size(elevels_ag_solv, 1), size(elevels_au_solv, 1), size(elevels_b3g_solv, 1), ...
        size(elevels_b3u_solv, 1)];
    elseif (l1_l3_same_parity == true)
        size_array(index_size, :) = [N, size(elevels_ag, 1), size(elevels_b1g, 1), ...
        size(elevels_b2g, 1), size(elevels_b3g, 1), size(elevels_au, 1), ...
        size(elevels_b2u, 1), size(elevels_b2u, 1), size(elevels_b3u, 1), ...
        size(elevels_ag_solv, 1), size(elevels_au_solv, 1), size(elevels_b2g_solv, 1), ...
        size(elevels_b2u_solv, 1)];
    elseif (l2_l3_same_parity == true)
        size_array(index_size, :) = [N, size(elevels_ag, 1), size(elevels_b1g, 1), ...
        size(elevels_b2g, 1), size(elevels_b3g, 1), size(elevels_au, 1), ...
        size(elevels_b1u, 1), size(elevels_b2u, 1), size(elevels_b3u, 1), ...
        size(elevels_ag_solv, 1), size(elevels_au_solv, 1), size(elevels_b1g_solv, 1), ...
        size(elevels_b1u_solv, 1)];
    end

    % Find proportion of solvable states
    if (all_same_parity == true) 
        solvable_prop(index_size, :) = [size(elevels_ag_solv, 1), size(elevels_au_solv, 1)] / (total_sites);
    elseif (l1_l2_same_parity == true)
        solvable_prop(index_size, :) = [size(elevels_ag_solv, 1), size(elevels_au_solv, 1), size(elevels_b3g_solv, 1), size(elevels_b3u_solv, 1)] / (total_sites);
    elseif (l1_l3_same_parity == true)
        solvable_prop(index_size, :) = [size(elevels_ag_solv, 1), size(elevels_au_solv, 1), size(elevels_b2g_solv, 1), size(elevels_b2u_solv, 1)] / (total_sites);
    elseif (l2_l3_same_parity == true)
        solvable_prop(index_size, :) = [size(elevels_ag_solv, 1), size(elevels_au_solv, 1), size(elevels_b1g_solv, 1), size(elevels_b1u_solv, 1)] / (total_sites);
    end

    %% Plot the level spacings histogram and show r values
    % Plot histograms
    figure(index_n)
    tiledlayout(2,4,'TileSpacing', 'tight','Padding','Tight')
    bin_factor = 5;
    
    nexttile
    histogram(spacings_ag, ceil(size(elevels_ag, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Ag')
    subtitle(['r = ' num2str(r_ag) '; total = ' num2str(size_array(index_size, 2))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b1g, ceil(size(elevels_b1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B1g')
    subtitle(['r = ' num2str(r_b1g) '; total = ' num2str(size_array(index_size, 3))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b2g, ceil(size(elevels_b2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B2g')
    subtitle(['r = ' num2str(r_b2g) '; total = ' num2str(size_array(index_size, 4))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b3g, ceil(size(elevels_b3g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B3g')
    subtitle(['r = ' num2str(r_b3g) '; total = ' num2str(size_array(index_size, 5))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_au, ceil(size(elevels_au, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('Au')
    subtitle(['r = ' num2str(r_au) '; total = ' num2str(size_array(index_size, 6))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b1u, ceil(size(elevels_b1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B1u')
    subtitle(['r = ' num2str(r_b1u) '; total = ' num2str(size_array(index_size, 7))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b2u, ceil(size(elevels_b2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B2u')
    subtitle(['r = ' num2str(r_b2u) '; total = ' num2str(size_array(index_size, 8))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_b3u, ceil(size(elevels_b3u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('B3u')
    subtitle(['r = ' num2str(r_b3u) '; total = ' num2str(size_array(index_size, 9))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    set(figure(index_n),'position',[0,100,2000,400])
    saveas(gcf, [folderName '/N=' num2str(N) '_elevel_hist.jpeg']);

    % Plot histogram of solvable wavefunctions
    if (all_same_parity == true) 
        figure(index_n + 1)
        tiledlayout(1,2, 'TileSpacing', 'tight','Padding','Tight')
        bin_factor = 5;
    
        nexttile
        histogram(spacings_ag_solv, ceil(size(elevels_ag_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('Ag Solvable')
        subtitle(['r = ' num2str(r_ag_solv) '; total = ' num2str(size(elevels_ag_solv, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    
        nexttile
        histogram(spacings_au_solv, ceil(size(elevels_au_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('Au Solvable')
        subtitle(['r = ' num2str(r_au_solv) '; total = ' num2str(size(elevels_au_solv, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    else
        figure(index_n + 1)
        tiledlayout(2, 2, 'TileSpacing', 'tight','Padding','Tight')
        bin_factor = 5;

        nexttile
        histogram(spacings_ag_solv, ceil(size(elevels_ag_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('Ag Solvable')
        subtitle(['r = ' num2str(r_ag_solv) '; total = ' num2str(size(elevels_ag_solv, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')
    
        nexttile
        histogram(spacings_au_solv, ceil(size(elevels_au_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title('Au Solvable')
        subtitle(['r = ' num2str(r_au_solv) '; total = ' num2str(size(elevels_au_solv, 1))])
        xlabel('Energy spacing')
        ylabel('Counts')

        if (l1_l2_same_parity == true)
        
            nexttile
            histogram(spacings_b3g_solv, ceil(size(elevels_b3g_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
            title('B3g Solvable')
            subtitle(['r = ' num2str(r_b3g_solv) '; total = ' num2str(size(elevels_b3g_solv, 1))])
            xlabel('Energy spacing')
            ylabel('Counts')
        
            nexttile
            histogram(spacings_b3u_solv, ceil(size(elevels_b3u_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
            title('B3u Solvable')
            subtitle(['r = ' num2str(r_b3u_solv) '; total = ' num2str(size(elevels_b3u_solv, 1))])
            xlabel('Energy spacing')
            ylabel('Counts')

        elseif (l1_l3_same_parity == true)
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

        elseif (l2_l3_same_parity == true)
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
    end

    set(figure(index_n + 1),'position',[0, 100, 700, 400])
    saveas(gcf, [folderName '/N=' num2str(N) '_elevel_hist_solv.jpeg']);
    
    % Increment indices
    index_n = index_n + 2;
    index_size = index_size + 1;
    index_r = index_r + 1;
end

%% Plot r as a function of size in each irrep
figure(index_n)
tiledlayout(2,4,'TileSpacing', 'tight','Padding','Tight')
ylow_unsolv = 0.45;
yhigh_unsolv = 0.7;
ylow_solv = 0.35;
yhigh_solv = 0.6;

nexttile
plot(N_nums, r_array(:, 2), 'Linestyle','-','Marker','.', [0 0.4470 0.7410])
title('Ag')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 3).', 'Linestyle','-','Marker','.', [0 0.4470 0.7410])
title('B1g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 4).', 'Linestyle','-','Marker','.', [0 0.4470 0.7410])
title('B2g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 5).', 'Linestyle','-','Marker','.', [0 0.4470 0.7410])
title('B3g')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 6).', 'Linestyle','-','Marker','.', [0 0.4470 0.7410])
title('Au')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 7).', 'Linestyle','-','Marker','.', [0 0.4470 0.7410])
title('B1u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 8).', 'Linestyle','-','Marker','.', [0 0.4470 0.7410])
title('B2u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

nexttile
plot(N_nums, r_array(:, 9).', 'Linestyle','-','Marker','.', [0 0.4470 0.7410])
title('B3u')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53')
ylim([ylow_unsolv, yhigh_unsolv])

set(figure(index_n),'position',[0,100,2000,400])
saveas(gcf, [folderName '/LSR_plot.jpeg']);

% Plot r values for solvable wavefunctions

if (all_same_parity == true) 
    figure(index_n + 1)
    tiledlayout(1, 2, 'TileSpacing', 'tight','Padding','Tight')

    nexttile
    plot(N_nums.', r_array(:, 10).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    title('Ag')
    xlabel('N')
    ylabel('r')
    hold on 
    yline(0.39, '--', 'r = 0.39')
    ylim([ylow_solv, yhigh_solv])

    nexttile
    plot(N_nums, r_array(:, 11).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    title('Au')
    xlabel('N')
    ylabel('r')
    hold on 
    yline(0.39, '--', 'r = 0.39')
    ylim([ylow_solv, yhigh_solv])

    set(figure(index_n + 1),'position',[0,100,1000,400])
    saveas(gcf, [folderName '/LSR_plot_solv.jpeg']);
 
else
    figure(index_n + 1)
    tiledlayout(2, 2,'TileSpacing', 'tight','Padding','Tight')

    nexttile
    plot(N_nums.', r_array(:, 10).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    title('Ag')
    xlabel('N')
    ylabel('r')
    hold on 
    yline(0.39, '--', 'r = 0.39')
    ylim([ylow_solv, yhigh_solv])

    nexttile
    plot(N_nums, r_array(:, 11).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
    title('Au')
    xlabel('N')
    ylabel('r')
    hold on 
    yline(0.39, '--', 'r = 0.39')
    ylim([ylow_solv, yhigh_solv])

    if (l1_l2_same_parity == true)
        nexttile
        plot(N_nums, r_array(:, 12).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
        title('B3g')
        xlabel('N')
        ylabel('r')
        hold on 
        yline(0.39, '--', 'r = 0.39')
        ylim([ylow_solv, yhigh_solv])
        
        nexttile
        plot(N_nums, r_array(:, 13).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
        title('B3u')
        xlabel('N')
        ylabel('r')
        hold on 
        yline(0.39, '--', 'r = 0.39')
        ylim([ylow_solv, yhigh_solv])
    elseif (l1_l3_same_parity == true)
        nexttile
        plot(N_nums, r_array(:, 12).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
        title('B1g')
        xlabel('N')
        ylabel('r')
        hold on 
        yline(0.39, '--', 'r = 0.39')
        ylim([ylow_solv, yhigh_solv])
        
        nexttile
        plot(N_nums, r_array(:, 13).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
        title('B1u')
        xlabel('N')
        ylabel('r')
        hold on 
        yline(0.39, '--', 'r = 0.39')
        ylim([ylow_solv, yhigh_solv])
    elseif (l2_l3_same_parity == true)
        nexttile
        plot(N_nums, r_array(:, 12).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
        title('B2g')
        xlabel('N')
        ylabel('r')
        hold on 
        yline(0.39, '--', 'r = 0.39')
        ylim([ylow_solv, yhigh_solv])
        
        nexttile
        plot(N_nums, r_array(:, 13).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
        title('B2u')
        xlabel('N')
        ylabel('r')
        hold on 
        yline(0.39, '--', 'r = 0.39')
        ylim([ylow_solv, yhigh_solv])
    end

    set(figure(index_n + 1),'position',[0,100,1000,400])
    saveas(gcf, [folderName '/LSR_plot_solv.jpeg']);
end

%% Plot r values together
figure(index_n + 2)
tiledlayout(1, 2,'TileSpacing', 'tight','Padding','Tight')
bin_factor = 5;
    
nexttile
plot(N_nums.', r_array(:, 2).', 'Linestyle','-','Marker','.')
hold on
plot(N_nums.', r_array(:, 3).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 4).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 5).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 6).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 7).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 8).', 'Linestyle','-','Marker','.')
plot(N_nums.', r_array(:, 9).', 'Linestyle','-','Marker','.')
legend('Ag','B1g','B2g','B3g','Au','B1u','B2u','B3u')
title('Unsolvable wavefunctions')
xlabel('N')
ylabel('r')
hold on 
yline(0.53, '--', 'r = 0.53', 'HandleVisibility','off')
ylim([ylow_unsolv, yhigh_unsolv])
hold off

nexttile
plot(N_nums, r_array(:, 10).', 'Linestyle','-','Marker','.')
hold on
plot(N_nums.', r_array(:, 11).', 'Linestyle','-','Marker','.')
if (all_same_parity == false)
    if (l1_l2_same_parity == true)
        plot(N_nums.', r_array(:, 12).', 'Linestyle','-','Marker','.')
        plot(N_nums.', r_array(:, 13).', 'Linestyle','-','Marker','.')
        legend('Ag', 'Au', 'B3g', 'B3u')
    elseif (l1_l3_same_parity == true)
        plot(N_nums.', r_array(:, 12).', 'Linestyle','-','Marker','.')
        plot(N_nums.', r_array(:, 13).', 'Linestyle','-','Marker','.')
        legend('Ag', 'Au', 'B1g', 'B1u')
    elseif (l2_l3_same_parity == true)
        plot(N_nums.', r_array(:, 12).', 'Linestyle','-','Marker','.')
        plot(N_nums.', r_array(:, 13).', 'Linestyle','-','Marker','.')
        legend('Ag', 'Au', 'B2g', 'B2u')
    end
else
    legend('Ag', 'Au')
end
title('Solvable wavefunctions')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39', 'HandleVisibility','off')
ylim([ylow_solv, yhigh_solv])
set(figure(index_n + 2),'position',[0,100,1000,400])
saveas(gcf, [folderName '/LSR_plot_combined.jpeg']);

%% Plot proportions of solvable states
figure(index_n + 3)
proportions = sum(solvable_prop, 2);
plot(N_nums, proportions, 'Linestyle','-','Marker','.')
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
save([folderName '/elevels_ag.mat'],'elevels_ag','-v7.3')
save([folderName '/elevels_au.mat'],'elevels_au','-v7.3')
save([folderName '/elevels_b1g.mat'],'elevels_b1g','-v7.3')
save([folderName '/elevels_b1u.mat'],'elevels_b1u','-v7.3')
save([folderName '/elevels_b2g.mat'],'elevels_b2g','-v7.3')
save([folderName '/elevels_b2u.mat'],'elevels_b2u','-v7.3')
save([folderName '/elevels_b3g.mat'],'elevels_b2g','-v7.3')
save([folderName '/elevels_b3u.mat'],'elevels_b2u','-v7.3')

save([folderName 'elevels_ag_solv_sorted.mat'],'elevels_ag_solv_sorted','-v7.3')

if (all_same_parity == true) 
    save([folderName '/elevels_au_solv_sorted.mat'],'elevels_au_solv_sorted','-v7.3')   
elseif (l1_l2_same_parity == true)
    save([folderName '/elevels_b3g_solv_sorted.mat'],'elevels_b3g_solv_sorted','-v7.3')
    save([folderName '/elevels_b3u_solv_sorted.mat'],'elevels_b3u_solv_sorted','-v7.3')
    save([folderName '/elevels_au_solv.mat'],'elevels_au_solv','-v7.3')
elseif (l1_l3_same_parity == true)
    save([folderName '/elevels_b2g_solv_sorted.mat'],'elevels_b2g_solv_sorted','-v7.3')
    save([folderName '/elevels_b2u_solv_sorted.mat'],'elevels_b2u_solv_sorted','-v7.3')
    save([folderName '/elevels_au_solv.mat'],'elevels_au_solv','-v7.3')
elseif (l2_l3_same_parity == true)
    save([folderName '/elevels_b1g_solv_sorted.mat'],'elevels_b1g_solv_sorted','-v7.3')
    save([folderName '/elevels_b1u_solv_sorted.mat'],'elevels_b1u_solv_sorted','-v7.3')
    save([folderName '/elevels_au_solv.mat'],'elevels_au_solv','-v7.3')
end

%% Functions
% Gives x coordinate of index within each face
function x = xfacecoord(index, N, l1, l2, l3, total_sites)
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
function y = yfacecoord(index, N, l1, l2, l3, total_sites) 
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
function i = index_from_coord(x, y, face, N, l1, l2, l3, total_sites)
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

% Given index on discretized cube returns new index after C2x symmetry
function i = c2x_index(index, N, l1, l2, l3, total_sites)
    if (index <= l1 * l3 * N^2) % in face 1
        new_face = 4;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index <= l1 * (l3 + l2) * N^2) % in face 2 
        new_face = 2;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l2 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index > total_sites - l1 * l2 * N^2) % in face 6
        new_face = 6;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l2 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % in face 3
        new_face = 5;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % in face 4
        new_face = 1;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    else  % in face 5
        new_face = 3;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    end
end

% Given index on discretized cube returns new index after inv symmetry
function i = inv_index(index, N, l1, l2, l3, total_sites)
    if (index <= l1 * l3 * N^2) % in face 1
        new_face = 4;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index <= l1 * (l3 + l2) * N^2) % in face 2  
        new_face = 6;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index > total_sites - l1 * l2 * N^2) % in face 6
        new_face = 2;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % in face 3
        new_face = 5;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % in face 4
        new_face = 1;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    else  % in face 5
        new_face = 3;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    end
end

% Given index on discretized cube returns new index after C2y symmetry
function i = c2y_index(index, N, l1, l2, l3, total_sites)
    if (index <= l1 * l3 * N^2) % in face 1
        new_face = 1;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index <= l1 * (l3 + l2) * N^2) % in face 2 
        new_face = 6;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l2 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index > total_sites - l1 * l2 * N^2) % in face 6
        new_face = 2;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l2 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % in face 3
        new_face = 5;
        new_x = l2 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % in face 4
        new_face = 4;
        new_x = l1 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    else  % in face 5
        new_face = 3;
        new_x = l2 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    end
end

% Given index on discretized cube returns new index after C2z symmetry
function i = c2z_index(index, N, l1, l2, l3, total_sites)
    if (index <= l1 * l3 * N^2) % in face 1
        new_face = 4;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index <= l1 * (l3 + l2) * N^2) % in face 2 
        new_face = 6;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (index > total_sites - l1 * l2 * N^2) % in face 6
        new_face = 2;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= l2 * N) % in face 3
        new_face = 3;
        new_x = l2 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    elseif (mod(index - l1 * (l3 + l2) * N^2 - 1, (2 * l2 + l1) * N) + 1 <= (l2 + l1) * N) % in face 4
        new_face = 1;
        new_x = xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    else  % in face 5
        new_face = 5;
        new_x = l2 * N + 1 - xfacecoord(index, N, l1, l2, l3, total_sites);
        new_y = l3 * N + 1 - yfacecoord(index, N, l1, l2, l3, total_sites);
        i = index_from_coord(new_x, new_y, new_face, N, l1, l2, l3, total_sites);
    end
end

% Functions to find indices of H matrix
function y = upper(x, N, l1, l2, l3, total_sites) 
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
        y = mod(x - (total_sites - l1 * l2 * N^2) - 1, l1 * N) + 1;
    else % x in non-top row upper arm or non-top row lower arm
        y = x + l1 * N;
    end
end

function y = lower(x, N, l1, l2, l3, total_sites) 
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

function y = right(x, N, l1, l2, l3, total_sites)
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

function y = left(x, N, l1, l2, l3, total_sites)
    if (x <= l3 * l1 * N^2 && (floor((x - 1) / (l1 * N)) == (x - 1) / (l1 * N))) % x in lower half of lower arm, leftmost column
        row_from_bottom = (x - 1) / (l1 * N) + 1;
        y = total_sites - l2 * l1 * N^2 - (row_from_bottom) * (2 * l2 + l1) * N + 1;
    elseif (x > l3 * l1 * N^2 && x <= (l3 + l2) * l1 * N^2 && floor((x - 1) / (l1 * N)) == (x - 1) / (l1 * N)) % x in upper half of lower arm, leftmost column
        row_from_bottom_2 = ceil(x / (l1 * N)) - l3 * N;
        y = (l3 + l2) * l1 * N^2 + row_from_bottom_2;
    elseif (x > (l3 + l2) * l1 * N^2 && x <= total_sites - l2 * l1 * N^2 && floor((x - 1 - (l3 + l2) * l1 * N^2)/((2 * l2 + l1) * N)) == (x - 1 - (l2 + l3) * l1 * N^2)/((2 * l2 + l1) * N)) % x in middle arm, leftmost column 
        row_from_top = ceil((total_sites - l2 * l1 * N^2 - x)/((2 * l2 + l1) * N));
        y = l1 * N * (row_from_top - 1) + 1;
    elseif (x > total_sites - l2 * l1 * N^2 && floor((x - (total_sites - l2 * l1 * N^2) - 1) / (l1 * N)) == (x - ((total_sites - l2 * l1 * N^2)) - 1) / (l1 * N))  % x in upper arm, leftmost column 
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