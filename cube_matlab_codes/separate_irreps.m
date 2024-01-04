N = 75;

elevels_a1 = zeros(6 * N^2, 5);
elevels_a2 = zeros(6 * N^2, 5);
elevels_e = zeros(6 * N^2, 5);
elevels_t1 = zeros(6 * N^2, 5);
elevels_t2 = zeros(6 * N^2, 5);

index_a1 = 1;
index_a2 = 1;
index_e = 1;
index_t1 = 1;
index_t2 = 1;

for i = 1:size(energy_levels, 1)
    if (energy_levels(i, 2) == 1) && (energy_levels(i, 3) == 1)
        elevels_a1(index_a1, :) = energy_levels(i, :);
        index_a1 = index_a1 + 1;
    elseif (energy_levels(i, 2) == 2) && (energy_levels(i, 4) == 2)
        elevels_a1(index_a1, :) = energy_levels(i, :);
        index_a1 = index_a1 + 1;
    elseif (energy_levels(i, 2) == 1) && (energy_levels(i, 3) == -1)
        elevels_a2(index_a2, :) = energy_levels(i, :);
        index_a2 = index_a2 + 1;
    elseif (energy_levels(i, 2) == 2) && (energy_levels(i, 4) == -1)
        elevels_e(index_e, :) = energy_levels(i, :);
        index_e = index_e + 1;
    elseif (energy_levels(i, 2) == 3) && (energy_levels(i, 3) == -1)
        elevels_t1(index_t1, :) = energy_levels(i, :);
        index_t1 = index_t1 + 1;
    elseif (energy_levels(i, 2) == 3) && (energy_levels(i, 3) == 1)
        elevels_t2(index_t2, :) = energy_levels(i, :);
        index_t2 = index_t2 + 1;
    else 
        fprintf('Error')
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

%% Find spacings of energy levels
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

%% Plot histograms
tiledlayout(1,5)
nexttile
histogram(spacings_a1, 200, 'FaceColor','black', 'EdgeColor','none')
title('A1')
xlabel('Energy spacing')
ylabel('Counts')

nexttile
histogram(spacings_a2, 50, 'FaceColor','black', 'EdgeColor','none')
title('A1')
xlabel('Energy spacing')
ylabel('Counts')

nexttile
histogram(spacings_e, 300, 'FaceColor','black', 'EdgeColor','none')
title('E')
xlabel('Energy spacing')
ylabel('Counts')

nexttile
histogram(spacings_t1, 400, 'FaceColor','black', 'EdgeColor','none')
title('T1')
xlabel('Energy spacing')
ylabel('Counts')

nexttile
histogram(spacings_t2, 400, 'FaceColor','black', 'EdgeColor','none')
title('T2')
xlabel('Energy spacing')
ylabel('Counts')

set(figure(1),'position',[0,300,1500,200])

