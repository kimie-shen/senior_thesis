% Plot histograms before and after separation

elevels_solv = elevels_ag_solv_sorted;
elevels_unsolv = elevels_ag;

elevels_all = cat(1, elevels_solv, elevels_unsolv);
elevels_all_sorted = sortrows(elevels_all);

% Compute spacings
spacings_solv = zeros(size(elevels_solv, 1) - 1, 1);
spacings_unsolv = zeros(size(elevels_unsolv, 1) - 1, 1);
spacings_all = zeros(size(elevels_all_sorted, 1) - 1, 1);

for i = 1:(size(elevels_solv, 1) - 1)
    spacings_solv(i) = abs(elevels_solv(i, 1) - elevels_solv(i + 1, 1));
end
for i = 1:(size(elevels_unsolv, 1) - 1)
    spacings_unsolv(i) = abs(elevels_unsolv(i, 1) - elevels_unsolv(i + 1, 1));
end
for i = 1:(size(elevels_all_sorted, 1) - 1)
    spacings_all(i) = abs(elevels_all_sorted(i, 1) - elevels_all_sorted(i + 1, 1));
end

% Compute LSR values
r_solv = LSR(spacings_solv(:, 1));
r_unsolv = LSR(spacings_unsolv(:, 1));
r_all = LSR(spacings_all(:, 1));

% Plot histogram
bin_factor = 5;

tiledlayout(1, 3, 'TileSpacing', 'tight','Padding','Tight')
nexttile
histogram(spacings_all(2:end), ceil(size(elevels_all_sorted, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title('Before Separation')
subtitle(['r = ' num2str(r_all) '; total = ' num2str(size(elevels_all_sorted, 1))])
xlabel('Energy spacing')
ylabel('Counts')

nexttile
histogram(spacings_solv(2:end), ceil(size(elevels_solv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title('After Separation: Solvable')
subtitle(['r = ' num2str(r_solv) '; total = ' num2str(size(elevels_solv, 1))])
xlabel('Energy spacing')
ylabel('Counts')

nexttile
histogram(spacings_unsolv(2:end), ceil(size(elevels_unsolv, 1)/bin_factor), 'FaceColor','black', 'EdgeColor','none')
title('After Separation: Unsolvable')
subtitle(['r = ' num2str(r_unsolv) '; total = ' num2str(size(elevels_unsolv, 1))])
xlabel('Energy spacing')
ylabel('Counts')

set(figure(1),'position',[0,100,1500,300])

function r = LSR(spacings)
    %Allocate variable for r_n values
    r_n = zeros(size(spacings, 1) - 1, 1);
    
    % Calculate r_n variables
    for i = 1:(size(spacings, 1) - 1)
        r_n(i) = min(spacings(i), spacings(i + 1)) ...
            / max(spacings(i), spacings(i + 1));
    end

    % Calculate r variable
    r = mean(r_n);
end
