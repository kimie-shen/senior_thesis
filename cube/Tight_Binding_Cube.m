clear
N = 75; % Number of sites per side
L = 1/N; % Spacing
H = zeros(6 * N^2, 6 * N^2); % Hamiltonian matrix
show_plot = false; % Plot (1) or no plot (0)

tic;
% Set matrix values
for i = 1:6*N^2
    H(i, i) = 4/L^2;
    H(i, upper(i, N)) = -1/L^2;
    H(i, lower(i, N)) = -1/L^2;
    H(i, right(i, N)) = -1/L^2;
    H(i, left(i, N)) = -1/L^2;
end

% Find eigenvalues
%[eigenvectors, eigenvalues] = eig(H);
toc


% Functions to find indices
function y = upper(x, N) 
    if (x > 2 * N^2 - N && x <= 2 * N^2) % x in top row of lower arm 
        y = 2 * N^2 + N + mod(x - 1, N) + 1;
    elseif (x > 2 * N^2 && x <= 5 * N^2 - 3 * N) % x in middle arm but not top row
        y = x + 3 * N;
    elseif (x > 5 * N^2 - 3 * N && x <= 5 * N^2 - 2 * N) % x top row of middle arms, left side
        dist_from_vertex = x - (5 * N^2 - 3 * N); % distance from leftmost vertex
        y = 6 * N^2 - dist_from_vertex * N + 1; 
    elseif (x > 5 * N^2 - 2 * N && x <= 5 * N^2 - 1 * N) % x in top row of middle arms, center side
        y = 5 * N^2 + mod(x - 1, N) + 1;
    elseif (x > 5 * N^2 - 1 * N && x <= 5 * N^2) % x in top row of middle arms, right side
        dist_from_fold = x - (5 * N^2 - 1 * N);
        y = 5 * N^2 + dist_from_fold * N;
    elseif (x > 5 * N^2 + N * (N - 1)) % x in top row of upper arm
        y = mod(x - 1, N) + 1;
    else % x in non-top row upper arm or non-top row lower arm
        y = x + N;
    end
end

function y = lower(x, N) 
    if (x > 0 && x <= N) % x in lower arm bottom row
        y = 6 * N^2 - N + mod(x - 1, N) + 1;
    elseif (x > 2 * N^2 && x <= 2 * N^2 + N) % x in middle arms, bottom row, left side
        dist_from_vertex = x - 2 * N^2;
        y = N^2 + (dist_from_vertex - 1) * N + 1;
    elseif (x > 2 * N^2 + N && x <= 2 * N^2 + 2 * N) % x in middle arms, bottom row, center
        y = x - 2 * N;
    elseif (x > 2 * N^2 + 2 * N && x <= 2 * N^2 + 3 * N) % x in middle arms, bottom row, right side
        dist_from_fold = x - (2 * N^2 + 2 * N);
        y = 2 * N^2 - (dist_from_fold - 1) * N ;
    elseif (x > 2 * N^2 + 3 * N && x <= 5 * N^2) % x in middle arms, above bottom row
        y = x - 3 * N;
    elseif (x > 5 * N^2 && x <= 5 * N^2 + N) % x in upper arm, bottom row
        y = x - 2 * N;
    else % x in upper arm, above bottom row or in lower arm above bottom row
        y = x - N;
    end
end

function y = right(x, N)
    if (x <= N^2 && floor(x / N) == x / N) % x in lower half of lower arm, rightmost column
        row_from_bottom = x/N;
        y = 5 * N^2 - (row_from_bottom - 1) * 3 * N;
    elseif (x > N^2 && x <= 2 * N^2 && floor(x / N) == x / N) % x in upper half of lower arm, rightmost column
        row_from_top = (2 * N^2 - x)/N + 1;
        y = 2 * N^2 + 2 * N + row_from_top;
    elseif (x > 2 * N^2 && x <= 5 * N^2 && floor((x - 2 * N^2)/(3 * N)) == (x - 2 * N^2)/(3 * N)) % x in middle arm, rightmost column 
        row_from_top_2 = (5 * N^2 - x)/ (3 * N) + 1;
        y = row_from_top_2 * N;
    elseif (x > 5 * N^2 && floor(x / N) == x / N)  % x in upper arm, rightmost column 
        row_from_bottom_2 = (x - 5 * N^2) / N;
        y = 5 * N^2 - N + row_from_bottom_2;
    else  
        y = x + 1;
    end
end

function y = left(x, N)
    if (x <= N^2 && floor((x - 1) / N) == (x - 1) / N) % x in lower half of lower arm, leftmost column
        row_from_bottom = (x - 1) / N + 1;
        y = 5 * N^2 - (row_from_bottom) * 3 * N + 1;
    elseif (x > N^2 && x <= 2 * N^2 && floor((x - 1) / N) == (x - 1) / N) % x in upper half of lower arm, rightmost column
        row_from_bottom_2 = ceil(x / N) - N;
        y = 2 * N^2 + row_from_bottom_2;
    elseif (x > 2 * N^2 && x <= 5 * N^2 && floor((x - 1 - 2 * N^2)/(3 * N)) == (x - 1 - 2 * N^2)/(3 * N)) % x in middle arm, leftmost column 
        row_from_top = ceil((5 * N^2 - x)/(3 * N));
        y = N * (row_from_top - 1) + 1;
    elseif (x > 5 * N^2 && floor((x - 1) / N) == (x - 1) / N)  % x in upper arm, rightmost column 
        row_from_top_2 = ceil((6 * N^2 - x) / N);
        y = 5 * N^2 - 3 * N + row_from_top_2;
    else  
        y = x - 1;
    end
end