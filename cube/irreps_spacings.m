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

%histogram(spacings,100)