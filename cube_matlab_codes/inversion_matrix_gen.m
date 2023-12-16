% Generate matrix for C3 symmetry about diagonal axis
N = 75;
inversion = zeros(6 * N^2, 6 * N^2);
for i = 1: (6 * N^2)
    inversion(i, inverted_index(i, N)) = 1;
end

% Gives x coordinate of index within each face
function x = xfacecoord(index, N)
    x = mod(index - 1, N) + 1;
end

% Gives y coordinate of index within each face
function y = yfacecoord(index, N) 
    if (index <= N^2) % in face 1
        y = ceil(index / N);
    elseif (index <= 2 * N^2) % in face 2 
        y = ceil((index - N^2) / N);
    elseif (index <= 5 * N^2) % in faces 3, 4, or 5
        y = ceil((index - 2 * N^2) / (3 * N));
    else % in face 6
        y = ceil((index - 5 * N^2) / N);
    end
end

% From face coordinate and face number give index
function i = index_from_coord(x, y, face, N)
    if (face == 1)
        i = x + N * (y - 1);
    elseif (face == 2)
        i = N^2 + x + N * (y - 1);
    elseif (face == 3)
        i = 2 * N^2 + x + 3 * N * (y - 1);
    elseif (face == 4)
        i = 2 * N^2 + N + x + 3 * N * (y - 1);
    elseif (face == 5)
        i = 2 * N^2 + 2 * N + x + 3 * N * (y - 1);
    else 
        i = 5 * N^2 + x + N * (y - 1);
    end
end

% From given index return index after rotation
function i = inverted_index(index, N)
    if (index <= N^2) % in face 1
        new_face = 4;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index <= 2 * N^2) % in face 2 
        new_face = 6;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (index > 5 * N^2) % in face 6
        new_face = 2;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 2 * N^2 - 1, 3 * N) + 1 <= N) % in face 3
        new_face = 5;
        new_x = xfacecoord(index, N);
        new_y = N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (mod(index - 2 * N^2 - 1, 3 * N) + 1 <= 2 * N) % in face 4
        new_face = 1;
        new_x = N + 1 - xfacecoord(index, N);
        new_y = yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    else  % in face 5
        new_face = 3;
        new_x = xfacecoord(index, N);
        new_y = N + 1 - yfacecoord(index, N);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end
