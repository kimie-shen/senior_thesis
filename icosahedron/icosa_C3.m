% Calculate r values for various N
clear;
close all
tolerance = 1e-8;  % Tolerance for degenerate energy levels
corner_potential = 0;
energy_cut_factor = 1;

% Find max N
upper_lim = 30000; % Max number of sites allowed in laptop memory
max_N = floor(sqrt(upper_lim) - 1);
fprintf(['Max N = ' num2str(max_N) '\n'])
N_nums = primes(max_N);
N_nums = N_nums(5:end);
%N_nums = [29, 31];

% Make directory
folderName = ['icosahedron_C5_plots_ecut=' num2str(energy_cut_factor)];
mkdir(folderName);

r_ag_array = zeros(size(N_nums, 2), 6);
r_au_array = zeros(size(N_nums, 2), 6);
size_array = zeros(size(N_nums, 1), 3);

index_n = 1;
index_nn = 1;
index_r = 1;
index_size = 1;

w = exp(-1i * 2 *pi / 5);
nn_vals = [0, 1, 2, 3, 4]; % 0, 1, 2, 3, 4

for n = 1:size(N_nums, 2)
    for j = 1:size(nn_vals, 2)
        N = N_nums(n);
        nn = nn_vals(j);
        total_sites = (N + 1)^2;
    
        %% Diagonalize H and sigma_d matrices
        tic;
        fprintf(['\n    N = ' num2str(N) ' n = ' num2str(nn) '\n'])
    
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
    
        % Edges
        for i = 1:total_sites
            if (on_right(i, N))
                H(i, right(i, N)) = H(i, right(i, N)) * w^nn;
            end
    
            if (on_left(i, N))
                 H(i, left(i, N)) = H(i, left(i, N)) * w^(-nn);
            end
        end
    
        % Calculate eigenvectors and eigenvalues
        [eigenvectors, eigenvalues] = eig(H);
        fprintf('evals and evecs done \n')
        toc
        
        %% Generate Symmetry Matrices
    
        % Generate matrix for inv
        inv = zeros(total_sites, total_sites);
        for i = 1: (total_sites)
            inv(i, inv_index(i, N)) = 1;
        end
    
        for i = 1:total_sites
            if (face_from_index(i, N) == 1)
                inv(i, inv_index(i, N)) = inv(i, inv_index(i, N)) * w^(2*nn);
            elseif (face_from_index(i, N) == 2)
                inv(i, inv_index(i, N)) = inv(i, inv_index(i, N)) * w^(2*nn);
            elseif (face_from_index(i, N) == 3)
                inv(i, inv_index(i, N)) = inv(i, inv_index(i, N)) * w^(3*nn);
            elseif (face_from_index(i, N) == 4)
                inv(i, inv_index(i, N)) = inv(i, inv_index(i, N)) * w^(3*nn);
            end
        end
        
        %% Calculate characters of symmetries
        tic;
        % Preallocate data array
        inv_evals = zeros(total_sites, 1);
    
        % Determine 1/3rd of max energy
        max_energy = max(diag(eigenvalues));
        max_energy_index = 0;
    
        for i = 1:total_sites
            if (eigenvalues(i, i) <= max_energy / energy_cut_factor)
                max_energy_index = max_energy_index + 1;
            end
        end
    
        for i = 1:max_energy_index
            inv_evals(i) = (conj(eigenvectors(:, i)).' * inv * eigenvectors(:, i));
        end
        fprintf('inv done \n')
        
        toc;
        
        %% Analyze Degeneracies
        % Reorder energies and evecs
        energy_levels = zeros(max_energy_index, 4);
        [energies, id] = sort(diag(eigenvalues));
        eigenvectors = eigenvectors(:, id);
        inv_evals_sorted = inv_evals(id);
        
        % Analyze for degeneracies
        energy = energies(1);
        degeneracy = 1;
        index = 1;
        trace_inv = inv_evals_sorted(1);
        
        for i = 2:max_energy_index
            if (abs(energies(i) - energy) < tolerance) % Next energy is degenerate
                degeneracy = degeneracy + 1;
                trace_inv = trace_inv + inv_evals_sorted(i);
            else % Next energy is new
                % Record stats for previous energy level
                energy_levels(index, 1) = energy;
                energy_levels(index, 2) = degeneracy;
                energy_levels(index, 3) = trace_inv;
                energy_levels(index, 4) = i;
        
                % Reset variables for new energy level
                energy = energies(i);
                degeneracy = 1;
                index = index + 1;
                trace_inv = inv_evals_sorted(i);
            end
        
            % Record energy level if we reach the end
            if (i == max_energy_index)
                energy_levels(index, 1) = energy;
                energy_levels(index, 2) = degeneracy;
                energy_levels(index, 3) = trace_inv;
                energy_levels(index, 4) = i;
            end
        end
        
        energy_levels = energy_levels(1:index, :, :);
        
        %% Separate the representation classes
        % Allocate space for each irrep
        elevels_ag = zeros(max_energy_index, 4);
        elevels_au = zeros(max_energy_index, 4);
        energy_levels_rounded = zeros(size(energy_levels, 1), 4);
        
        index_ag = 1;
        index_au = 1;
        
        % Round the energy_levels characters    
        energy_levels_rounded(:, 3) = round(energy_levels(:, 3));
        energy_levels_rounded(:, 1) = energy_levels(:, 1);
        energy_levels_rounded(:, 2) = energy_levels(:, 2);
        energy_levels_rounded(:, 4) = energy_levels(:, 4);
    
    
        % Fill in elevel info in each irrep
        for i = 1:size(energy_levels_rounded, 1)
            trace_e = energy_levels_rounded(i, 2);
            trace_inv = energy_levels_rounded(i, 3);
        
            traces = [trace_e, trace_inv];
        
            if (isequal(traces, [1, 1])) % Ag
                elevels_ag(index_ag, :) = energy_levels_rounded(i, :);
                index_ag = index_ag + 1;
            elseif (isequal(traces, [1, -1])) % Au
                elevels_au(index_au, :) = energy_levels_rounded(i, :);
                index_au = index_au + 1;
            else 
                fprintf([num2str(energy_levels_rounded(i, 1)) ' ' num2str(energy_levels_rounded(i, 2)) ' ' ...
                    num2str(energy_levels_rounded(i, 3)) ' ' num2str(energy_levels_rounded(i, 4)) '\n'])
            end
        end
        
        % Remove extra rows of zeros
        if (index_ag > 1)
            elevels_ag = elevels_ag(1:(index_ag - 1), :);
        end 
    
        if (index_au > 1)
            elevels_au = elevels_au(1:(index_au - 1), :);
        end 

        %% Find energy level spacings
        % Make vector of energy level spacings
        if ((j == 1) && (n == 1))
            spacings_ag = zeros(size(elevels_ag, 1) - 1, 1);
            spacings_au = zeros(size(elevels_au, 1) - 1, 1);
        end
        
        for i = 1:(size(elevels_ag, 1) - 1)
            spacings_ag(i) = abs(elevels_ag(i, 1) - elevels_ag(i+1, 1));
        end
               
        for i = 1:(size(elevels_au, 1) - 1)
            spacings_au(i) = abs(elevels_au(i, 1) - elevels_au(i+1, 1));
        end
    
        % Compute level spacing ratios
        r_ag = LSR(spacings_ag);
        r_au = LSR(spacings_au);
    
        % Fill in r array values
        r_ag_array(index_r, 1) = N;
        r_au_array(index_r, 1) = N;
        r_ag_array(index_r, j + 1) = r_ag;
        r_au_array(index_r, j + 1) = r_au;
       
       
        % Fill in size of irreps
        size_array(index_size, :) = [N, size(elevels_ag, 1), size(elevels_au, 1)];

        %% Plot the level spacings histogram and show r values
        % Plot histograms
        figure(index_n)
        if (j == 1)
            tiledlayout(2, 5,'TileSpacing', 'tight','Padding','Tight')
        end
        bin_factor = 5;
        
        nexttile(j)
        histogram(spacings_ag, ceil(size(elevels_ag, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title(['Ag, n = ' num2str(nn)])
        subtitle(['r = ' num2str(r_ag) '; total = ' num2str(size_array(index_size, 2))])
        xlabel('Energy spacing')
        ylabel('Counts')
    
        nexttile(j + 5)
        histogram(spacings_au, ceil(size(elevels_au, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
        title(['Au, n = ' num2str(nn)])
        subtitle(['r = ' num2str(r_au) '; total = ' num2str(size_array(index_size, 3))])
        xlabel('Energy spacing')
        ylabel('Counts')
    
        set(figure(index_n),'position',[0,100,3000,400])
        saveas(gcf, [folderName '/N=' num2str(N) '_elevel_hist.jpeg']);
        
        % Increment indices
        index_size = index_size + 1;
    end

    index_n = index_n + 1;
    index_r = index_r + 1;
end




%% Plot r as a function of size in each irrep
figure(index_n + 1)
tiledlayout(2, 5, 'TileSpacing', 'tight', 'Padding', 'Tight')
ylow = 0.3;
yhigh = 0.6;

nexttile
plot(N_nums.', r_ag_array(:, 2).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Ag, n = 0')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_au_array(:, 2).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Au, n = 0')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_ag_array(:, 3).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Ag, n = 1')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_au_array(:, 3).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Au, n = 1')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_ag_array(:, 4).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Ag, n = 2')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums.', r_au_array(:, 4).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Au, n = 2')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_ag_array(:, 5).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Ag, n = 3')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_au_array(:, 5).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Au, n = 3')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_ag_array(:, 6).', 'Linestyle', '-', 'Marker', '.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Ag, n = 4')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

nexttile
plot(N_nums, r_au_array(:, 6).', 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
title('Au, n = 5')
xlabel('N')
ylabel('r')
hold on 
yline(0.39, '--', 'r = 0.39')
ylim([ylow, yhigh])

set(figure(index_n + 1),'position',[0,100,3000,400])
saveas(gcf, [folderName '/LSR_plot.jpeg']);

%% Plot proportions of solvable states
%figure(index_n + 1)
%plot(N_nums, sum(solvable_prop, 2), 'Linestyle','-','Marker','.', 'MarkerEdgeColor', [0 0.4470 0.7410])
%title('Proportion of states which are solvable')
%xlabel('N')
%yline(1/12, '--', '1/12')
%ylim([0.0, 0.1])
%ylabel('Proportion')
%saveas(gcf, [folderName '/solv_prop.jpeg']);

%% Save variables
save([folderName '/r_ag_array.mat'],'r_ag_array','-v7.3')
save([folderName '/r_au_array.mat'],'r_au_array','-v7.3')
save([folderName '/size_array.mat'],'size_array','-v7.3')
save([folderName '/elevels_ag.mat'],'elevels_ag','-v7.3')
save([folderName '/elevels_au.mat'],'elevels_au','-v7.3')


%% Functions


% Gives width of row within face CHECKED
function w = row_width(index, N) 
    face = face_from_index(index, N);
    y = yfacecoord(index, N);
    if ((face == 1) || (face == 3))
        w = 1 + 2 * (y - 1);
    else
        w = N - 2 * (y - 1);
    end
end

% Gives y coordinate of index within each face CHECKED
function y = yfacecoord(index, N) 
    sites_per_face = ((N + 1) / 2)^2;
    face = face_from_index(index, N);
    if (face == 1 || face == 3)
        index = mod(index - 1, sites_per_face) + 1;
        y = ceil(sqrt(index));
    else
        index = mod(index - 1, sites_per_face) + 1;
        y = (N + 1) / 2 + 1 - ceil(sqrt(sites_per_face + 1 - index));
    end
end

% Gives x coordinate of index within each face CHECKED
function x = xfacecoord(index, N) 
    face = face_from_index(index, N);
    y_coord = yfacecoord(index, N);
    sites_per_face = ((N + 1) / 2)^2;
    row_w = row_width(index, N);

    if (face == 1 || face == 3) % bottom row
        index = mod(index - 1, sites_per_face) + 1;
        x = index - (y_coord - 1)^2; 
    elseif (face == 2)
        total_sites = 2 * sites_per_face;
        x = row_w - (total_sites + 1 - index - ((N + 1) / 2 - y_coord)^2) + 1;
    elseif (face == 4)
        total_sites = 4 * sites_per_face;
        x = row_w - (total_sites + 1 - index - ((N + 1) / 2 - y_coord)^2) + 1;
    end
end

% Gives face containing given index CHECKED
function f = face_from_index(index, N) 
    sites_per_face = ((N + 1) / 2)^2;
    if (index <= 1 * sites_per_face)
        f = 1;
    elseif (index <= 2 * sites_per_face)
        f = 2;
    elseif (index <= 3 * sites_per_face)
        f = 3;
    else
        f = 4;
    end
end

% From face coordinate and face number give index CHECKED
function i = index_from_coord(x, y, face, N)
    sites_per_face = ((N + 1) / 2)^2;

    if (face == 1 || face == 3)
        i = sites_per_face * (face - 1) + (y - 1)^2 + x;
    else
        total_sites = face * sites_per_face;
        i = total_sites - ((N + 1) / 2 + 1 - y)^2 + x;
    end
end

% Given index on discretized octahedron returns new index after inv CHECKED
function i = inv_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    if (face == 1) 
        new_face = 4;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 3;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 2;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 1;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Functions to find indices of H matrix 
% CHECKED
function y = upper(index, N) 
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face == 1 && y_coord < (N + 1) / 2) % in face 1 but not top row
        y = index_from_coord(x_coord + 1, y_coord + 1, 1, N);
    elseif (face == 1 && y_coord == (N + 1) / 2) % in face 1 top row
        y = index_from_coord(x_coord, 1, 2, N);
    elseif (face == 3 && y_coord < (N + 1) / 2) % in face 3 but not top row
        y = index_from_coord(x_coord + 1, y_coord + 1, 3, N);
    elseif (face == 3 && y_coord == (N + 1) / 2) % in face 3 top row
        y = index_from_coord(x_coord, 1, 4, N);
    else % faces 2, 4
        y = index_from_coord(x_coord - 1, y_coord + 1, face, N);
    end
end

% CHECKED
function y = lower(index, N) 
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face == 1 || face == 3)
        y = index_from_coord(x_coord - 1, y_coord - 1, face, N);
    elseif (face == 2 && y_coord == 1) % face 2 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 1, N);
    elseif (face == 4 && y_coord == 1) % face 4 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 3, N);
    else
        y = index_from_coord(x_coord + 1, y_coord - 1, face, N);
    end
end

% CHECKED
function y = right(index, N) 
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);

    if (face == 1 && x_coord == row_w) % face 1 right side
        y = index_from_coord(1, y_coord, 1, N);
    elseif (face == 2 && x_coord == row_w) % face 2 right side
        y = index_from_coord(1, y_coord, 3, N);
    elseif (face == 3 && x_coord == row_w) % face 3 right side
        y = index_from_coord(1, y_coord, 2, N);
    elseif (face == 4 && x_coord == row_w) % face 4 right side
        y = index_from_coord(1, y_coord, 4, N);
    else
        y = index + 1; 
    end
end

% CHECKED
function y = left(index, N) 
    face = face_from_index(index, N);
    row_w = row_width(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face == 1 && x_coord == 1) % face 1 left side
        y = index_from_coord(row_w, y_coord, 1, N);
    elseif (face == 2 && x_coord == 1) % face 2 left side
        y = index_from_coord(N + 1 - row_w, y_coord, 3, N);
    elseif (face == 3 && x_coord == 1) % face 3 left side
        y = index_from_coord(N + 1 - row_w, y_coord, 2, N);
    elseif (face == 4 && x_coord == 1) % face 4 left side
        y = index_from_coord(row_w, y_coord, 4, N);
    else
        y = index - 1;
    end  
end

% Indicate whether nearest vertical neighbor is above or below
function answer = has_upper(index, N)
    y = yfacecoord(index, N);
    face = face_from_index(index, N);
    sites_per_face = ((N + 1) / 2)^2;

    if (face == 1 || face == 3) % in lower row 
        index = mod(index - 1, sites_per_face) + 1; 
        if (mod(y, 2) - mod(index, 2) == 0)
            answer = true;
        else 
            answer = false;
        end
    else
        index = index - (face - 1) * sites_per_face;
        if (mod(y, 2) - mod(index, 2) == 0)
            answer = false;
        else 
            answer = true;
        end
    end
end

% Indicate whether index is on right edge CHECKED
function answer = on_right(index, N)
    x = xfacecoord(index, N);
    row_w = row_width(index, N);
    face = face_from_index(index, N);

    if (face == 2)
        answer = false;
    elseif (x == row_w)
        answer = true;
    else
        answer = false;
    end
end

% Indicate whether index is on left edge CHECKED
function answer = on_left(index, N)
    x = xfacecoord(index, N);
    face = face_from_index(index, N);

    if (face == 3)
        answer = false;
    elseif (x == 1)
        answer = true;
    else
        answer = false;
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