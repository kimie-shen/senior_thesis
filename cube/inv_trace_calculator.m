% Preallocate column vector of i characters
inv_trace = zeros(size(energy_levels, 1), 1);

index = 1; 

for i = 1:(size(energy_levels, 1))
    % Sum over degenerate energy levels
    next_index = index + energy_levels(i, 2) - 1;
    inv_trace(i) = sum(inv_exp(index:next_index)); 

    % Update index
    index = index + energy_levels(i, 2);
end


