%% Separate energy_levels by irreps in O_h representation
N = 75;

% Round energy_levels variable
energy_levels(:, 2:6) = round(energy_levels(:, 2:6));

elevels_a1g = zeros(6 * N^2, 7);
elevels_a2g = zeros(6 * N^2, 7);
elevels_eg = zeros(6 * N^2, 7);
elevels_t1g = zeros(6 * N^2, 7);
elevels_t2g = zeros(6 * N^2, 7);
elevels_a1u = zeros(6 * N^2, 7);
elevels_a2u = zeros(6 * N^2, 7);
elevels_eu = zeros(6 * N^2, 7);
elevels_t1u = zeros(6 * N^2, 7);
elevels_t2u = zeros(6 * N^2, 7);

index_a1g = 1;
index_a2g = 1;
index_eg = 1;
index_t1g = 1;
index_t2g = 1;
index_a1u = 1;
index_a2u = 1;
index_eu = 1;
index_t1u = 1;
index_t2u = 1;

for i = 1:size(energy_levels, 1)
    trace_e = energy_levels(i, 2);
    trace_sigma_d = energy_levels(i, 3);
    trace_c3 = energy_levels(i, 4);
    trace_inv = energy_levels(i, 5);
    trace_c2 = energy_levels(i, 6);

    traces = [trace_e, trace_sigma_d, trace_c3, trace_inv, trace_c2];

    if (isequal(traces, [1, 1, 1, 1, 1])) % A1g
        elevels_a1g(index_a1g, :) = energy_levels(i, :);
        index_a1g = index_a1g + 1;
    elseif (isequal(traces, [2, 0, 2, 0, 2])) % Accidental degeneracy A1g + A1u
        elevels_a1g(index_a1g, :) = energy_levels(i, :);
        index_a1g = index_a1g + 1;
        elevels_a1u(index_a1u, :) = energy_levels(i, :); 
        index_a1u = index_a1u + 1;
    elseif (isequal(traces, [2, 0, 2, 0, -2])) % Accidental degeneracy A2g + A2u
        elevels_a2g(index_a2g, :) = energy_levels(i, :);
        index_a2g = index_a2g + 1;
        elevels_a2u(index_a2u, :) = energy_levels(i, :);
        index_a2u = index_a2u + 1;
    elseif (isequal(traces, [3, 1, 3, 1, -1])) % Accidental degeneracy A1g + A2u + A2g
        elevels_a1g(index_a1g, :) = energy_levels(i, :);
        index_a1g = index_a1g + 1;
        elevels_a2u(index_a2u, :) = energy_levels(i, :);
        index_a2u = index_a2u + 1;
        elevels_a2g(index_a2g, :) = energy_levels(i, :);
        index_a2g = index_a2g + 1;
    elseif (isequal(traces, [3, -1, 3, 1, 1])) % Accidental degeneracy A1g + A2g + A1u
        elevels_a1g(index_a1g, :) = energy_levels(i, :);
        index_a1g = index_a1g + 1;
        elevels_a2g(index_a2g, :) = energy_levels(i, :);
        index_a2g = index_a2g + 1;
        elevels_a1u(index_a1u, :) = energy_levels(i, :);
        index_a1u = index_a1u + 1;
    elseif (isequal(traces, [1, -1, 1, 1, -1])) % A2g
        elevels_a2g(index_a2g, :) = energy_levels(i, :);
        index_a2g = index_a2g + 1;
    elseif (isequal(traces, [2, 0, -1, 2, 0])) % Eg
        elevels_eg(index_eg, :) = energy_levels(i, :);
        index_eg = index_eg + 1;
    elseif (isequal(traces, [3, -1, 0, 3, -1])) % T1g
        elevels_t1g(index_t1g, :) = energy_levels(i, :);
        index_t1g = index_t1g + 1;
    elseif (isequal(traces, [3, 1, 0, 3, 1])) %T2g
        elevels_t2g(index_t2g, :) = energy_levels(i, :);
        index_t2g = index_t2g + 1;
    elseif (isequal(traces, [1, -1, 1, -1, 1])) %A1u
        elevels_a1u(index_a1u, :) = energy_levels(i, :);
        index_a1u = index_a1u + 1;
    elseif (isequal(traces, [1, 1, 1, -1, -1])) %A2u
        elevels_a2u(index_a2u, :) = energy_levels(i, :);
        index_a2u = index_a2u + 1;
    elseif (isequal(traces, [2, 0, -1, -2, 0])) %Eu
        elevels_eu(index_eu, :) = energy_levels(i, :);
        index_eu = index_eu + 1;
    elseif (isequal(traces, [3, 1, 0, -3, -1])) %T1u
        elevels_t1u(index_t1u, :) = energy_levels(i, :);
        index_t1u = index_t1u + 1;
    elseif (isequal(traces, [3, -1, 0, -3, +1])) %T2u
        elevels_t2u(index_t2u, :) = energy_levels(i, :);
        index_t2u = index_t2u + 1;
    else 
        fprintf([num2str(energy_levels(i, 1)) ' ' num2str(energy_levels(i, 2)) ' ' ...
            num2str(energy_levels(i, 3)) ' ' num2str(energy_levels(i, 4)) ' ' ...
            num2str(energy_levels(i, 5)) ' ' num2str(energy_levels(i, 6)) ' ' ...
            num2str(energy_levels(i, 7)) '\n'])
    end
end

% Remove extra rows of zeros
if (index_a1g > 1)
    elevels_a1g = elevels_a1g(1:(index_a1g - 1), :);
end 

if (index_a2g > 1)
    elevels_a2g = elevels_a2g(1:(index_a2g - 1), :);
end

if (index_eg > 1)
    elevels_eg = elevels_eg(1:(index_eg - 1), :);
end

if (index_t1g > 1)
    elevels_t1g = elevels_t1g(1:(index_t1g - 1), :);
end

if (index_t2g > 1)
    elevels_t2g = elevels_t2g(1:(index_t2g - 1), :);
end

if (index_a1u > 1)
    elevels_a1u = elevels_a1u(1:(index_a1u - 1), :);
end 

if (index_a2u > 1)
    elevels_a2u = elevels_a2u(1:(index_a2u - 1), :);
end

if (index_eu > 1)
    elevels_eu = elevels_eu(1:(index_eu - 1), :);
end

if (index_t1u > 1)
    elevels_t1u = elevels_t1u(1:(index_t1u - 1), :);
end

if (index_t2u > 1)
    elevels_t2u = elevels_t2u(1:(index_t2u - 1), :);
end

%% Find energy level spacings
% Make vector of energy level spacings
spacings_a1g = zeros(size(elevels_a1g, 1) - 1, 1);
spacings_a2g = zeros(size(elevels_a2g, 1) - 1, 1);
spacings_eg = zeros(size(elevels_eg, 1) - 1, 1);
spacings_t1g = zeros(size(elevels_t1g, 1) - 1, 1);
spacings_t2g = zeros(size(elevels_t2g, 1) - 1, 1);
spacings_a1u = zeros(size(elevels_a1u, 1) - 1, 1);
spacings_a2u = zeros(size(elevels_a2u, 1) - 1, 1);
spacings_eu = zeros(size(elevels_eu, 1) - 1, 1);
spacings_t1u = zeros(size(elevels_t1u, 1) - 1, 1);
spacings_t2u = zeros(size(elevels_t2u, 1) - 1, 1);

for i = 1:(size(elevels_a1g, 1) - 1)
    spacings_a1g(i) = abs(elevels_a1g(i, 1) - elevels_a1g(i+1, 1));
end

for i = 1:(size(elevels_a2g, 1) - 1)
    spacings_a2g(i) = abs(elevels_a2g(i, 1) - elevels_a2g(i+1, 1));
end

for i = 1:(size(elevels_eg, 1) - 1)
    spacings_eg(i) = abs(elevels_eg(i, 1) - elevels_eg(i+1, 1));
end

for i = 1:(size(elevels_t1g, 1) - 1)
    spacings_t1g(i) = abs(elevels_t1g(i, 1) - elevels_t1g(i+1, 1));
end

for i = 1:(size(elevels_t2g, 1) - 1)
    spacings_t2g(i) = abs(elevels_t2g(i, 1) - elevels_t2g(i+1, 1));
end

for i = 1:(size(elevels_a1u, 1) - 1)
    spacings_a1u(i) = abs(elevels_a1u(i, 1) - elevels_a1u(i+1, 1));
end

for i = 1:(size(elevels_a2u, 1) - 1)
    spacings_a2u(i) = abs(elevels_a2u(i, 1) - elevels_a2u(i+1, 1));
end

for i = 1:(size(elevels_eu, 1) - 1)
    spacings_eu(i) = abs(elevels_eu(i, 1) - elevels_eu(i+1, 1));
end

for i = 1:(size(elevels_t1u, 1) - 1)
    spacings_t1u(i) = abs(elevels_t1u(i, 1) - elevels_t1u(i+1, 1));
end

for i = 1:(size(elevels_t2u, 1) - 1)
    spacings_t2u(i) = abs(elevels_t2u(i, 1) - elevels_t2u(i+1, 1));
end


% Compute level spacing ratios
r_a1g = LSR(spacings_a1g);
r_a1u = LSR(spacings_a1u);
r_a2g = LSR(spacings_a2g);
r_a2u = LSR(spacings_a2u);
r_eg = LSR(spacings_eg);
r_eu = LSR(spacings_eu);
r_t1g = LSR(spacings_t1g);
r_t1u = LSR(spacings_t1u);
r_t2g = LSR(spacings_t2g);
r_t2u = LSR(spacings_t2u);

% Plot histograms
figure(1)
tiledlayout(2,5)
bin_factor = 5;
    
nexttile
histogram(spacings_a1g, ceil(size(elevels_a1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title(['A1g: r = ' num2str(r_a1g)])
xlabel('Energy spacing')
ylabel('Counts')
    
nexttile
histogram(spacings_a2g, ceil(size(elevels_a2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title(['A2g: r = ' num2str(r_a2g)])
xlabel('Energy spacing')
ylabel('Counts')
    
nexttile
histogram(spacings_eg, ceil(size(elevels_eg, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title(['Eg: r = ' num2str(r_eg)])
xlabel('Energy spacing')
ylabel('Counts')
    
nexttile
histogram(spacings_t1g, ceil(size(elevels_t1g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title(['T1g: r = ' num2str(r_t1g)])
xlabel('Energy spacing')
ylabel('Counts')
    
nexttile
histogram(spacings_t2g, ceil(size(elevels_t2g, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title(['T2g: r = ' num2str(r_t2g)])
xlabel('Energy spacing')
ylabel('Counts')
    
nexttile
histogram(spacings_a1u, ceil(size(elevels_a1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title(['A1u: r = ' num2str(r_a1u)])
xlabel('Energy spacing')
ylabel('Counts')
    
nexttile
histogram(spacings_a2u, ceil(size(elevels_a2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title(['A2u: r = ' num2str(r_a2u)])
xlabel('Energy spacing')
ylabel('Counts')
    
nexttile
histogram(spacings_eu, ceil(size(elevels_eu, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title(['Eu: r = ' num2str(r_eu)])
xlabel('Energy spacing')
ylabel('Counts')

nexttile
histogram(spacings_t1u, ceil(size(elevels_t1u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title(['T1u: r = ' num2str(r_t1u)])
xlabel('Energy spacing')
ylabel('Counts')
    
nexttile
histogram(spacings_t2u, ceil(size(elevels_t2u, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title(['T2u: r = ' num2str(r_t2u)])
xlabel('Energy spacing')
ylabel('Counts')
    
set(figure(1),'position',[0,100,1500,400])


function r = LSR(spacings)
    r_n = zeros(size(spacings, 1) - 1, 1);
    
    for i = 1:(size(spacings, 1) - 1)
        r_n(i) = min(spacings(i), spacings(i + 1)) ...
            / max(spacings(i), spacings(i + 1));
    end

    r = mean(r_n);
end