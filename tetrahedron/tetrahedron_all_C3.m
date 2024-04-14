% Calculate r values for various N
clear;
close all
tolerance = 1e-8;  % Tolerance for degenerate energy levels
corner_potential = 0;
e_cut_factor = 3;

% Find max N
upper_lim = 30000; % Max number of sites allowed in laptop memory
max_N = floor(sqrt(upper_lim) - 1);
fprintf(['Max N = ' num2str(max_N) '\n'])
N_nums = primes(max_N);
N_nums = N_nums(4:end);
%N_nums = [31];

r_array = zeros(size(N_nums, 1), 6);
size_array = zeros(size(N_nums, 1), 6);
solvable_prop = zeros(size(N_nums, 1), 2);

% Make directory
folderName = ['tetrahedron_plots_ecut=' num2str(e_cut_factor)];
mkdir(folderName);

index_n = 1;
index_r = 1;
index_size = 1;

for n = 1:size(N_nums, 2)
    N = N_nums(n);
    total_sites = (N + 1)^2;

    %% Diagonalize H and sigma_d matrices
    tic;
    fprintf(['\n    N = ' num2str(N) '\n'])

    % Generate Hamiltonian matrix
    L = 1/N; % Spacing
    onSiteHopping = 4 / L^2;
    neighborHopping = 4 / (3 * L^2);
    H = zeros((N + 1)^2, (N + 1)^2); % Hamiltonian matrix
    
    for i = 1:(N + 1)^2
        %fprintf([num2str(i) ' '])
        H(i, i) = onSiteHopping;
        %fprintf([num2str(right(i, N)) ' '])
        H(i, right(i, N)) = -neighborHopping;
        %fprintf([num2str(left(i, N)) ' '])
        H(i, left(i, N)) = -neighborHopping;
        if (has_upper(i, N))
            %fprintf('upper ')
        end
        if (has_upper(i, N))
            %fprintf([num2str(upper(i, N)) ' \n'])
            H(i, upper(i, N)) = -neighborHopping;
        else
            %fprintf([num2str(lower(i, N)) ' \n'])
            H(i, lower(i, N)) = -neighborHopping;
        end
    end

    %% Add corner potential
    faces = [1, 3, 4];
    for i = 1:3
        face = faces(i);

        lower_left_corner = index_from_coord(1, 1, face, N); 
        lower_right_corner = index_from_coord(N, 1, face, N); 
        top_corner = index_from_coord(1, (N + 1) / 2, face, N);

        %fprintf([num2str(lower_left_corner) ' ' num2str(lower_right_corner) ' ' num2str(top_corner) '\n'])

        H(lower_left_corner, lower_left_corner) = H(lower_left_corner, lower_left_corner) + corner_potential; 
        H(lower_right_corner, lower_right_corner) = H(lower_right_corner, lower_right_corner) + corner_potential;
        H(top_corner, top_corner) = H(top_corner, top_corner) + corner_potential;
    end

    % Face 2
    bottom_corner = index_from_coord(1, 1, 2, N);
    top_right_corner = index_from_coord(N, (N + 1) / 2, 2, N);
    top_left_corner = index_from_coord(1, (N + 1) / 2, 2, N);

    H(bottom_corner, bottom_corner) = H(bottom_corner, bottom_corner) + corner_potential; 
    H(top_right_corner, top_right_corner) = H(top_right_corner, top_right_corner) + corner_potential;
    H(top_left_corner, top_left_corner) = H(top_left_corner, top_left_corner) + corner_potential;

    %fprintf([num2str(bottom_corner) ' ' num2str(top_right_corner) ' ' num2str(top_left_corner) '\n'])


    % Generate matrix for C3
    C3 = zeros((N + 1)^2, (N + 1)^2);
    
    for i = 1:((N + 1)^2)
        %fprintf([num2str(i) ' ' num2str(c3_index(i, N)) '\n'])
        C3(i, c3_index(i, N)) = 1;
    end

    % Generate matrix for C3_inv
    C3_inv = zeros((N + 1)^2, (N + 1)^2);
    
    for i = 1:((N + 1)^2)
        %fprintf([num2str(i) ' ' num2str(c3_index(i, N)) '\n'])
        C3_inv(c3_index(i, N), i) = 1;
    end
    
    % Generate matrix for C2
    C2 = zeros((N + 1)^2, (N + 1)^2);
    
    for i = 1:((N + 1)^2)
        %fprintf([num2str(i) ' ' num2str(c2_index(i, N)) '\n'])
        C2(i, c2_index(i, N)) = 1;
    end

    % Generate matrix for sigma_d symmetry 
    sigma_d = zeros((N + 1)^2, (N + 1)^2);
    for i = 1: (N + 1)^2
        %fprintf([num2str(i) ' ' num2str(sigma_d_index(i, N)) '\n'])
        sigma_d(i, sigma_d_index(i, N)) = 1;
    end
    
    % Calculate eigenvectors and eigenvalues
    [eigenvectors, eigenvalues] = eig(H + ((1 + sqrt(5))/2000) * sigma_d + pi/1000 * (C3 + C3_inv));
    fprintf('evals and evecs done \n')
    toc
    
    
    %% Calculate characters of symmetries
    tic;
    % Preallocate data array
    sigma2_evals = zeros((N + 1)^2, 1);
    c3_evals = zeros((N + 1)^2, 1);
    c3_c3_inv_evals = zeros((N + 1)^2, 1);
    c2_evals = zeros((N + 1)^2, 1);
    
    
    for i = 1:((N + 1)^2)
        sigma2_evals(i) = (eigenvectors(:, i).' * sigma_d * eigenvectors(:, i));
        %eigenvalues(i, i) = eigenvalues(i, i) - ((1 + sqrt(5))/2000) * sigma2_evals(i);
    end
    fprintf('sigma_d done \n')
    
    for i = 1:((N + 1)^2)
        c3_evals(i) = (eigenvectors(:, i).' * C3 * eigenvectors(:, i));
    end
    fprintf('c3 done \n')

    for i = 1:((N + 1)^2)
        c3_c3_inv_evals(i) = (eigenvectors(:, i).' * (C3 + C3_inv)/2 * eigenvectors(:, i));
    end
    fprintf('c3 + c3_inv done \n')

    for i = 1:((N + 1)^2)
        c2_evals(i) = (eigenvectors(:, i).' * C2 * eigenvectors(:, i));
    end
    fprintf('c2 done \n')

    toc;

     % Determine 1/3rd of max energy
    max_energy = eigenvalues(total_sites, total_sites);
    max_energy_index = 0;

    for i = 1:total_sites
        if (eigenvalues(i, i) < max_energy / e_cut_factor)
            max_energy_index = max_energy_index + 1;
        end
    end
    
    %% Analyze Degeneracies
    % Reorder energies and evecs
    energy_levels = zeros(max_energy_index, 7);
    [energies, id] = sort(diag(eigenvalues));
    eigenvectors = eigenvectors(:, id);
    sigma2_evals_sorted = sigma2_evals(id);
    c3_evals_sorted = c3_evals(id);
    c3_c3_inv_evals_sorted = c3_c3_inv_evals(id);
    c2_evals_sorted = c2_evals(id);
    
    % Analyze for degeneracies
    energy = energies(1);
    degeneracy = 1;
    index = 1;
    trace_sigma2 = sigma2_evals_sorted(1);
    trace_c3 = c3_evals_sorted(1);
    trace_c3_c3_inv = c3_c3_inv_evals_sorted(1);
    trace_c2 = c2_evals_sorted(1);
    
    for i = 2:max_energy_index
        if (abs(energies(i) - energy) < tolerance) % Next energy is degenerate
            degeneracy = degeneracy + 1;
            trace_sigma2 = trace_sigma2 + sigma2_evals_sorted(i);
            trace_c3 = trace_c3 + c3_evals_sorted(i);
            trace_c3_c3_inv = trace_c3_c3_inv + c3_c3_inv_evals_sorted(i);
            trace_c2 = trace_c2 + c2_evals_sorted(i);
        else % Next energy is new
            % Record stats for previous energy level
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_sigma2;
            energy_levels(index, 4) = trace_c3;
            energy_levels(index, 5) = trace_c2;
            energy_levels(index, 6) = trace_c3_c3_inv;
            energy_levels(index, 7) = i - 1;
    
            % Reset variables for new energy level
            energy = energies(i);
            degeneracy = 1;
            index = index + 1;
            trace_sigma2 = sigma2_evals_sorted(i);
            trace_c3 = c3_evals_sorted(i);
            trace_c3_c3_inv = c3_c3_inv_evals_sorted(i);
            trace_c2 = c2_evals_sorted(i);
        end
    
        % Record energy level if we reach the end
        if (i == max_energy_index)
            energy_levels(index, 1) = energy;
            energy_levels(index, 2) = degeneracy;
            energy_levels(index, 3) = trace_sigma2;
            energy_levels(index, 4) = trace_c3;
            energy_levels(index, 5) = trace_c2;
            energy_levels(index, 6) = trace_c3_c3_inv;
            energy_levels(index, 7) = i;
        end
    end
    
    energy_levels = energy_levels(1:index, :);
    
    %% Separate the representation classes
    % Allocate space for each irrep
    elevels_a1 = zeros(max_energy_index, 6);
    elevels_a2 = zeros(max_energy_index, 6);
    elevels_e = zeros(max_energy_index, 6);
    elevels_t1 = zeros(max_energy_index, 6);
    elevels_t2 = zeros(max_energy_index, 6);
    
    index_a1 = 1;
    index_a2 = 1;
    index_e = 1;
    index_t1 = 1;
    index_t2 = 1;
    
    % Round the energy_levels characters
    energy_levels_rounded = zeros(size(energy_levels, 1), 6);
    energy_levels_rounded(:, 2:5) = round(energy_levels(:, 2:5));
    energy_levels_rounded(:, 1) = energy_levels(:, 1);
    energy_levels_rounded(:, 6) = energy_levels(:, 6);
    
    % Fill in elevel info in each irrep
    for i = 1:size(energy_levels_rounded, 1)
        trace_e = energy_levels_rounded(i, 2);
        trace_sigma_d = energy_levels_rounded(i, 3);
        trace_c3 = energy_levels_rounded(i, 4);
        trace_c2 = energy_levels_rounded(i, 5);
    
        traces = [trace_e, trace_sigma_d, trace_c3, trace_c2];
    
        if (isequal(traces, [1, 1, 1, 1])) % A1
            elevels_a1(index_a1, :) = energy_levels_rounded(i, :);
            index_a1 = index_a1 + 1;
        elseif (isequal(traces, [1, -1, 1, 1])) % A2
            elevels_a2(index_a2, :) = energy_levels_rounded(i, :);
            index_a2 = index_a2 + 1;
        elseif (isequal(traces, [2, 0, -1, 2])) % E
            elevels_e(index_e, :) = energy_levels_rounded(i, :);
            index_e = index_e + 1;
        elseif (isequal(traces, [3, -1, 0, -1])) % T1
            elevels_t1(index_t1, :) = energy_levels_rounded(i, :);
            index_t1 = index_t1 + 1;
        elseif (isequal(traces, [3, 1, 0, -1])) %T2
            elevels_t2(index_t2, :) = energy_levels_rounded(i, :);
            index_t2 = index_t2 + 1;
        elseif (isequal(traces, [3, -1, 0, 3])) % Accidental degeneracy E + A2
            elevels_e(index_e, :) = energy_levels_rounded(i, :);
            index_e = index_e + 1;
            elevels_a2(index_a2, :) = energy_levels_rounded(i, :); 
            index_a2 = index_a2 + 1;
        else 
            fprintf([num2str(energy_levels_rounded(i, 1)) ' ' num2str(energy_levels_rounded(i, 2)) ' ' ...
                num2str(energy_levels_rounded(i, 3)) ' ' num2str(energy_levels_rounded(i, 4)) ' ' ...
                num2str(energy_levels_rounded(i, 5)) '\n'])
        end
    end
    
    % Remove extra rows of zeros
    if (index_a1 > 1)
        elevels_a1 = elevels_a1(1:(index_a1 - 1), :);
    end 
    
    if (index_a2 > 1)
        elevels_a2 = elevels_a2(1:(index_a2 - 1), :);
    end
    
    if (index_e > 1)
        elevels_e = elevels_e(1:(index_e - 1), :);
    end
    
    if (index_t1 > 1)
        elevels_t1 = elevels_t1(1:(index_t1 - 1), :);
    end
    
    if (index_t2 > 1)
        elevels_t2 = elevels_t2(1:(index_t2 - 1), :);
    end
    
  
    
    %% Find energy level spacings
    % Make vector of energy level spacings
    spacings_a1 = zeros(size(elevels_a1, 1) - 1, 1);
    spacings_a2 = zeros(size(elevels_a2, 1) - 1, 1);
    spacings_e = zeros(size(elevels_e, 1) - 1, 1);
    spacings_t1 = zeros(size(elevels_t1, 1) - 1, 1);
    spacings_t2 = zeros(size(elevels_t2, 1) - 1, 1);
    
    for i = 1:(size(elevels_a1, 1) - 1)
        spacings_a1(i) = abs(elevels_a1(i, 1) - elevels_a1(i+1, 1));
    end
    
    for i = 1:(size(elevels_a2, 1) - 1)
        spacings_a2(i) = abs(elevels_a2(i, 1) - elevels_a2(i+1, 1));
    end
    
    for i = 1:(size(elevels_e, 1) - 1)
        spacings_e(i) = abs(elevels_e(i, 1) - elevels_e(i+1, 1));
    end
    
    for i = 1:(size(elevels_t1, 1) - 1)
        spacings_t1(i) = abs(elevels_t1(i, 1) - elevels_t1(i+1, 1));
    end
    
    for i = 1:(size(elevels_t2, 1) - 1)
        spacings_t2(i) = abs(elevels_t2(i, 1) - elevels_t2(i+1, 1));
    end
    
    % Compute level spacing ratios
    r_a1 = LSR(spacings_a1);
    r_a2 = LSR(spacings_a2);
    r_e = LSR(spacings_e);
    r_t1 = LSR(spacings_t1);
    r_t2 = LSR(spacings_t2);

    % Fill in r array values
    r_array(index_r, :) = [N, r_a1, r_a2, r_e, r_t1, r_t2];
   
    % Fill in size of irreps
    size_array(index_size, :) = [N, size(elevels_a1, 1), size(elevels_a2, 1), size(elevels_e, 1), ...
        size(elevels_t1, 1), size(elevels_t2, 1)];

    % Find proportion of solvable states
    solvable_prop(index_size, :) = [size(elevels_a1, 1), size(elevels_a2, 1)] / ((N + 1)^2);

    %% Plot the level spacings histogram and show r values
    % Plot histograms
    figure(index_n)
    tiledlayout(1, 5,'TileSpacing', 'tight','Padding','Tight')
    bin_factor = 5;
    
    nexttile
    histogram(spacings_a1, ceil(size(elevels_a1, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A1')
    subtitle(['r = ' num2str(r_a1) '; total = ' num2str(size_array(index_size, 2))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_a2, ceil(size(elevels_a2, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('A2')
    subtitle(['r = ' num2str(r_a2) '; total = ' num2str(size_array(index_size, 3))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_e, ceil(size(elevels_e, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('E')
    subtitle(['r = ' num2str(r_e) '; total = ' num2str(size_array(index_size, 4))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t1, ceil(size(elevels_t1, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('T1')
    subtitle(['r = ' num2str(r_t1) '; total = ' num2str(size_array(index_size, 5))])
    xlabel('Energy spacing')
    ylabel('Counts')
    
    nexttile
    histogram(spacings_t2, ceil(size(elevels_t2, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
    title('T2')
    subtitle(['r = ' num2str(r_t2) '; total = ' num2str(size_array(index_size, 6))])
    xlabel('Energy spacing')
    ylabel('Counts')

    set(figure(index_n),'position',[0,100,1500,200])
    saveas(gcf, [folderName '/N=' num2str(N) '_elevel_hist.jpeg']);
    
    % Increment indices
    index_n = index_n + 1;
    index_size = index_size + 1;
    index_r = index_r + 1;
end

%% Plot r as a function of size in each irrep
figure(index_n)
tiledlayout(1, 5, 'TileSpacing', 'tight', 'Padding', 'Tight')
ylow = 0.3;
yhigh = 0.7;

nexttile
plot(N_nums.', r_array(:, 2).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A1')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 3).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('A2')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 4).', 'Linestyle', '-','Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('E')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 5).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('T1')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_array(:, 6).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('T2')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

set(figure(index_n),'position',[0,100,1500,200])
saveas(gcf, [folderName '/LSR_plot.jpeg']);

% Save variables
save([folderName '/r_array.mat'],'r_array','-v7.3')
save([folderName '/size_array.mat'],'size_array','-v7.3')
save([folderName '/size_array.mat'],'size_array','-v7.3')
save([folderName '/elevels_a1.mat'],'elevels_a1','-v7.3')
save([folderName '/elevels_a2.mat'],'elevels_a2','-v7.3')
save([folderName '/elevels_e.mat'],'elevels_e','-v7.3')
save([folderName '/elevels_t1.mat'],'elevels_t1','-v7.3')
save([folderName '/elevels_t2.mat'],'elevels_t2','-v7.3')

%% Functions

% Gives width of row within face
function w = row_width(index, N) 
    face = face_from_index(index, N);
    y = yfacecoord(index, N);
    if (face == 2)
        w = 1 + 2 * (y - 1);
    else
        w = N - 2 * (y - 1);
    end
end

% Gives lnegth of right slanting column within each face
function r = right_col(x_coord, y_coord, N)
    matrix = zeros(N, (N + 1) / 2);
    for i = 1:N
        for j = 1:(N + 1) / 2
            matrix(i, j) = (1 + 2 * (j - 1)) + 2 * floor(i / 2);
        end
    end
    r = matrix(x_coord, y_coord);
end

% Gives length of right slanting column within each face
function r = left_col(x_coord, y_coord, N)
    matrix = zeros(N, (N + 1) / 2);
    for j = 1:(N + 1) / 2
        for i = 1:N
            matrix(i, j) = N - 2 * floor((i - 1) / 2);
        end
    end
    r = matrix(x_coord, y_coord);
end

% Gives face as a function of index
function f = face_from_index(index, N)
    sites_per_face = ((N + 1) / 2)^2;
    total_sites = (N + 1)^2;

    if (index > 3 * sites_per_face)
        f = 4;
    else 
        y = yfacecoord(index, N);
        outer_triangle_width = N - 2 * (y - 1);
        sites_per_row = 2 * N + 1 - 2 * (y - 1);

        x_all = index - (total_sites - (N + 2 - y)^2);

        %fprintf(['index ' num2str(index) '\n'])
        %fprintf(['y ' num2str(y) '\n'])
        %fprintf(['outer_triangle_width ' num2str(outer_triangle_width) '\n'])
        %fprintf(['sites_per_row ' num2str(sites_per_row) '\n'])
        %fprintf(['x_all ' num2str(x_all) '\n \n'])

        if (x_all <= outer_triangle_width)
            f = 1;
        elseif (x_all > sites_per_row - outer_triangle_width)
            f = 3;
        else
            f= 2;
        end
    end
end

% Gives x coordinate of index within each face
function x = xfacecoord(index, N)
    rows_all = N + 2 - ceil(sqrt((N + 1)^2 + 1 - index));
    total_sites = (N + 1)^2;
    x_all = index - (total_sites - (N + 2 - rows_all)^2);
    face = face_from_index(index, N);

    if (face == 4)
        x = x_all;
    elseif (face == 1)
        x = x_all;
    elseif (face == 2)
        x = x_all - (N - 2 * (rows_all - 1));
    elseif (face == 3)
        x = x_all - (N + 1);
    end
end

% Gives y coordinate of index within each face
function y = yfacecoord(index, N) 
    rows_per_triangle = (N + 1) / 2;
    rows_all = N + 2 - ceil(sqrt((N + 1)^2 + 1 - index));
    y = mod(rows_all - 1, rows_per_triangle) + 1; 
end

% From face coordinate and face number give index
function i = index_from_coord(x, y, face, N)
    total_sites = (N + 1)^2;
    total_rows = N + 1;
    rows_above = total_rows - y;
    outer_triangle_width = N - 2 * (y - 1);

    %fprintf([' \nx ' num2str(x) '\n'])
    %fprintf(['y ' num2str(y) '\n'])
    %fprintf(['face ' num2str(face) '\n'])
    %fprintf(['rows_above ' num2str(rows_above) '\n'])
    %fprintf(['rows_above ' num2str((N + 1) / 2 - y) '\n'])
    %fprintf(['outer_triangle_width ' num2str(outer_triangle_width) '\n '])

    if (face == 1)
        i = total_sites - (rows_above + 1)^2 + x;
    elseif (face == 2)
        i = total_sites - (rows_above + 1)^2 + outer_triangle_width + x;
    elseif (face == 3)
        i = total_sites - (rows_above + 1)^2 + (N + 1) + x;
    elseif (face == 4)
        i = total_sites - (((N + 1) / 2) + 1  - y)^2 + x;
    end
end

% Given index on discretized cube returns new index after sigma_d symmetry
function i = sigma_d_index(index, N)
    face = face_from_index(index, N);
    row_w = row_width(index, N);
    if (face == 1) 
        new_face = 3;
        new_x = row_w + 1- xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 2;
        new_x = row_w + 1- xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 1;
        new_x = row_w + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 4;
        new_x = row_w + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized cube returns new index after C3 symmetry
function i = c3_index(index, N)
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    face = face_from_index(index, N);
    if (face == 1) % in face 1
        new_face = 3;
        new_x = left_col(x_coord, y_coord, N) + 1 - (2 * y_coord);
        new_y = ceil((xfacecoord(index, N))/ 2);
        %fprintf(['left_col = ' num2str(left_col(x_coord, y_coord, N)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        if (has_upper(index, N))
            i = index_from_coord(new_x, new_y, new_face, N);
        else
            i = index_from_coord(new_x, new_y, new_face, N) + 1;
        end
    elseif (face == 2) % in face 2 
        new_face = 2;
        new_y = (N + 1) / 2 + 1 - y_coord + floor(x_coord / 2);
        if (has_upper(index, N))
            new_x =  N + 1 - 2 * (y_coord - 1) - 1;
        else
            new_x =  N + 1 - 2 * (y_coord - 1);
        end
        %fprintf(['right_col = ' num2str(1 + 2 * (new_y - 1)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N) + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3)
        new_face = 4;
        new_x = left_col(x_coord, y_coord, N) + 1 - (2 * y_coord);
        new_y = ceil((xfacecoord(index, N))/ 2);
        %fprintf(['left_col = ' num2str(left_col(x_coord, y_coord, N)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        if (has_upper(index, N))
            i = index_from_coord(new_x, new_y, new_face, N);
        else
            i = index_from_coord(new_x, new_y, new_face, N) + 1;
        end
    elseif (face == 4) % in face 3
        new_face = 1;
        new_x = left_col(x_coord, y_coord, N) + 1 - (2 * y_coord);
        new_y = ceil((xfacecoord(index, N))/ 2);
        %fprintf(['left_col = ' num2str(left_col(x_coord, y_coord, N)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        if (has_upper(index, N))
            i = index_from_coord(new_x, new_y, new_face, N);
        else
            i = index_from_coord(new_x, new_y, new_face, N) + 1;
        end
    end
end

% Given index on discretized cube returns new index after C2 symmetry
function i = c2_index(index, N)
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    face = face_from_index(index, N);
    if (face == 1) % in face 1
        new_face = 3;
        new_x = x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) % in face 2 
        new_face = 4;
        row_w = row_width(index, N);
        new_y = (N + 1) / 2 + 1 - y_coord;
        new_x =  row_w + 1 - x_coord;
        %fprintf(['right_col = ' num2str(1 + 2 * (new_y - 1)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3)
        new_face = 1;
        new_x = x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) % in face 3
        new_face = 2;
        row_w = row_width(index, N);
        new_y = (N + 1) / 2 + 1 - y_coord;
        new_x =  row_w + 1 - x_coord;
        %fprintf(['left_col = ' num2str(left_col(x_coord, y_coord, N)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Functions to find indices of H matrix
function y = upper(x, N) 
    face = face_from_index(x, N);
    y_coord = yfacecoord(x, N);
    if (face == 4)
        y_coord = y_coord + (N + 1) / 2;
    end
    y = x + 2 * (N + 1 - y_coord);
end

function y = lower(x, N) 
    face = face_from_index(x, N);

    if (yfacecoord(x, N) == 1 && face == 1)
        y = index_from_coord(N + 1 - x, 1, 3, N);
        %fprintf(['y=' num2str(y)])
    elseif (yfacecoord(x, N) == 1 && face == 3)
        x_coord = xfacecoord(x, N);
        y = index_from_coord(N + 1 - x_coord, 1, 1, N);
    elseif (face == 4)
        y = x - 2 * ((N + 1) / 2 + 1 - yfacecoord(x, N));
    else
        y = x - 2 * (N + 2 - yfacecoord(x, N));
    end
end

function y = right(x, N)
    face = face_from_index(x, N);
    x_coord = xfacecoord(x, N);
    y_coord = yfacecoord(x, N);
    row_w = row_width(x, N);

    if (x_coord == row_w && face == 4) % top half, right edge
        y = index_from_coord(N + 1 - x_coord, (N + 1) / 2 + 1 - y_coord, 3, N);
    elseif (x_coord == row_w && face == 3) % bottom half, right edge
        %fprintf(['x_coord = ' num2str(x_coord) '\n'])
        %fprintf(['y_coord = ' num2str(y_coord) '\n'])
        %fprintf(['row_w = ' num2str(row_w) '\n'])
        %fprintf(['new_x_coord = ' num2str(N + 1 - x_coord) '\n'])
        %fprintf(['new_y_coord = ' num2str((N + 1) / 2 + 1 - y_coord) '\n'])
        y = index_from_coord(N + 1 - x_coord, (N + 1) / 2 + 1 - y_coord, 4, N);
    else
        y = x + 1;
    end
end

function y = left(x, N)
    face = face_from_index(x, N);
    row_w = row_width(x, N);
    x_coord = xfacecoord(x, N);
    y_coord = yfacecoord(x, N);

    if (x_coord == 1 && face == 4) % top half, left edge
        y = index_from_coord(1, ((N + 1) / 2 + 1 - y_coord), 1, N);
    elseif (x_coord == 1 && face == 1) % bottom half, left edge
        y = index_from_coord(1, ((N + 1) / 2 + 1 - y_coord), 4, N);
    else
        y = x - 1;
    end
end

% Indicate whether nearest vertical neighbor is above or below
function answer = has_upper(x, N)
    y = yfacecoord(x, N);
    face = face_from_index(x, N);

    if (face == 4)
        y = y + (N + 1) / 2;
    end

    if (mod(y, 2) - mod(x, 2) == 0)
        answer = false;
    else 
        answer = true;
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