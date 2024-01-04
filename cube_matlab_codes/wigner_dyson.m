
% Make vector of energy level spacings
spacings = zeros(6 * N^2 - 1, 1);
for i=1:(6*N^2 - 1)
    spacings(i) = (abs(e(i)) - abs(e(i+1)));
end

%histogram(spacings,100)

% Bin the energy level spacings
max_spacing = max(spacings);
num_bins = 100000;
bin_width = max_spacing/num_bins;
fprintf(num2str(bin_width))

binned_spacing = zeros(num_bins, 1);

for i = 1:(6*N^2 - 1)
    bin_num = ceil(spacings(i)/bin_width);

    if (bin_num > 0)
        binned_spacing(bin_num) = binned_spacing(bin_num) + 1;
    else
        binned_spacing(1) = binned_spacing(1) + 1;
    end
    
end

% Create vector for x-axis
energy_spacings = bin_width:bin_width:(num_bins * bin_width);

% Remove degenerate spacings
energy_spacings = energy_spacings(2:end);
binned_spacing = binned_spacing(2:end);

% Plot histogram
bar(energy_spacings, binned_spacing)
xlabel('Energy Spacing')
ylabel('Count')

