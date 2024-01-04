% Make vector of energy level spacings
spacings = zeros(6 * N^2 - 1);
for i=1:(6*N^2 - 1)
    spacings(i) = (abs(e(i)) - abs(e(i+1)));
end

%histogram(spacings,100)