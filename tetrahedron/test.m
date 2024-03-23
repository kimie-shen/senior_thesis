N = 5;
total = 2 * (N + 1)^2;

%x = [1, 1, 2, 3, 1, 2, 3, 4, 5];
%y = [1, 2, 2, 2, 3, 3, 3, 3, 3];

x = [1, 2, 3, 4, 5, 1, 2, 3, 1];
y = [1, 1, 1, 1, 1, 2, 2, 2, 3];

%index = [1, 2, 4, 5, 7, 9, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56, 58, 60, 62, 65, 67, 70];
index = [3, 6, 8, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 63, 64, 66, 68, 69, 71, 72];

for i = 1:total
    fprintf([num2str(i) ' '  mat2str(inv_index(i, N)) '\n'])
end

% Given index on discretized octahedron returns new index after C3 symmetry
% CHECKED
function i = inv_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    if (face == 1) 
        new_face = 7;
        new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 5;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = y_coord - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 1;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = (N + 1) / 2 - ceil(x_coord / 2) + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 3;
        new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 5) 
        new_face = 2;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 6) 
        new_face = 8;
        new_x = x_coord;
        new_y = (N + 1) / 2 - y_coord + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 7) 
        new_face = 6;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = (N + 1) / 2 - ceil(x_coord / 2) + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 8) 
        new_face = 4;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end


function i = c3_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    if (face == 1) 
        new_face = 6;
        new_x = x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 1;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 + 1 - y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 5;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 4;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 + 1 - y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 5) 
        new_face = 8;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 6) 
        new_face = 2;
        new_x = 1 + 2 * (y_coord - 1) - x_coord + 1;
        new_y = (N + 1) / 2 + 1 - ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 7) 
        new_face = 7;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = (N + 1) / 2 + 1 - y_coord - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 8) 
        new_face = 3;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized cube returns new index after C2 symmetry
function i = c2_index(index, N)
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    face = face_from_index(index, N);
    row_w = row_width(index, N);
    if (face == 1) % in face 1
        new_face = 8;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) % in face 2 
        new_face = 7;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        %fprintf(['right_col = ' num2str(1 + 2 * (new_y - 1)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3)
        new_face = 6;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) % in face 4
        new_face = 5;
        row_w = row_width(index, N);
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        %fprintf(['left_col = ' num2str(left_col(x_coord, y_coord, N)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 5) % in face 5
        new_face = 4;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 6) % in face 6 
        new_face = 3;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        %fprintf(['right_col = ' num2str(1 + 2 * (new_y - 1)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 7)
        new_face = 2;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 8) % in face 8
        new_face = 1;
        new_x = row_w + 1 - x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        %fprintf(['left_col = ' num2str(left_col(x_coord, y_coord, N)) ' index = ' num2str(index) ' new_x = ' num2str(new_x) ' new_y = ' num2str(new_y) '\n'])
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized cube returns new index after sigma_d symmetry
function i = sigma_d_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    if (face == 1) 
        new_face = 3;
        new_x = 2 * (y_coord - 1) + 1 - x_coord + 1;
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 7;
        new_x = (N + 1) - (2 * (y_coord - 1) + 1) + mod(x_coord - 1, 2);
        new_y = y_coord - floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 1;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 5;
        new_x = 2 * (y_coord - 1) + 1 - (x_coord - 1);
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 5) 
        new_face = 4;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 6) 
        new_face = 8;
        new_x = 2 * (y_coord - 1) + 1 - (x_coord - 1);
        new_y = ceil(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 7) 
        new_face = 2;
        new_x = N + 1 - (2 * (y_coord - 1) + 1) - x_coord + 1;
        new_y = (N + 1) / 2 - ceil(x_coord / 2) + 1;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 8) 
        new_face = 6;
        new_x = 2 * (y_coord - 1) + 1 + mod(x_coord - 1, 2);
        new_y = y_coord + floor(x_coord / 2);
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end


function y = lower(index, N) 
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);

    if (face == 1)
        y = index_from_coord(x_coord - 1, y_coord - 1, 1, N);
    elseif (face == 3 && y_coord == 1) % face 3 bottom row
        new_y = (x_coord - 1) / 2 + 1;
        y = index_from_coord(1, new_y, 1, N);
    elseif (face == 5 && y_coord == 1) % face 5 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 1, N);
    elseif (face == 7 && y_coord == 1) % face 7 bottom row
        new_y = (N + 1) / 2 + 1 - ((x_coord - 1) / 2 + 1); 
        y = index_from_coord(N + 1 - x_coord, new_y, 1, N);
    elseif (face == 8 && y_coord == 1) % face 8 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 4, N);
    elseif (face == 8)
        y = index_from_coord(x_coord, y_coord - 1, 8, N) + 1;
    else
        y = index - 3 * (N + 1) - 1;
    end
end

% Indicate whether nearest vertical neighbor is above or below
function answer = has_upper(index, N)
    y = yfacecoord(index, N);
    face = face_from_index(index, N);

    if (face == 1 || face == 8)
        if (mod(y, 2) - mod(index, 2) == 0)
            answer = true;
        else 
            answer = false;
        end
    elseif (mod(index, 2) == 0)
        answer = true;
    else
        answer = false;
    end
end


function y = right(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);

    if (face == 1 && x_coord == row_w) % face 1 right side
        new_x = N + 1 - (1 + 2 * (y_coord - 1));
        y = index_from_coord(new_x, 1, 7, N);
    elseif (face == 8 && x_coord == row_w) % face 8 right side
        new_x = 1 + 2 * (y_coord - 1);
        y = index_from_coord(new_x, (N + 1) / 2, 6, N);
    elseif (face == 7 && x_coord == row_w) % face 7 right side
        y = index_from_coord(1, y_coord, 2, N);
    else
        y = index + 1; 
    end
end

% Gives x coordinate of index within each face CHECKED
function x = xfacecoord(index, N) 
    face = face_from_index(index, N);
    y_coord = yfacecoord(index, N);
    sites_per_face = ((N + 1) / 2)^2;
    total_sites = 2 * (N + 1)^2;
    row_w = row_width(index, N);

    if (face == 1)
        x = index - (y_coord - 1)^2; 
    elseif (face == 8)
        x = row_w - (total_sites + 1 - index - ((N + 1) / 2 - y_coord)^2) + 1;
    elseif (face == 2) 
        x = index - sites_per_face - (y_coord - 1) * 3 * (N + 1);
    elseif (face == 4)
        x = index - sites_per_face - (y_coord - 1) * 3 * (N + 1) - (N + 1);
    elseif (face == 6)
        x = index - sites_per_face - (y_coord - 1) * 3 * (N + 1) - 2 * (N + 1);
    elseif (face == 3)
        x = index - sites_per_face - (y_coord - 1) * 3 * (N + 1) - (N + 1 - row_w);
    elseif (face == 5)
        x = index - sites_per_face - (y_coord - 1) * 3 * (N + 1) - (N + 1) - (N + 1 - row_w);
    elseif (face == 7)
        x = index - sites_per_face - (y_coord - 1) * 3 * (N + 1) - 2 * (N + 1) - (N + 1 - row_w);
    end
end



% Functions to find indices of H matrix
function y = upper(index, N) 
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);

    if (face == 1 && y_coord < (N + 1) / 2) % in face 1 but not top row
        y = index_from_coord(x_coord + 1, y_coord + 1, 1, N);
    elseif (face == 1 && y_coord == (N + 1) / 2) % in face 1 top row
        y = index_from_coord(x_coord, 1, 5, N);
    elseif (face == 8)
        y = index_from_coord(x_coord - 1, y_coord + 1, 8, N);
    elseif (face == 2 && y_coord == (N + 1) / 2) % face 2 top row
        new_y = (N + 1) / 2 - (x_coord - 1) / 2;
        y = index_from_coord(1, new_y, 8, N);
    elseif (face == 6 && y_coord == (N + 1) / 2) % face 6 top row
        new_y = (x_coord - 1) / 2 + 1;
        y = index_from_coord(row_w + 1 - x_coord, new_y, 8, N);
    elseif (face == 4 && y_coord == (N + 1) / 2) % face 4 top row
        y = index + 2 * (N + 1);
    elseif (mod(face, 2) == 0) % face = 2, 4, or 6
        y = index_from_coord(x_coord + 1, y_coord + 1, face, N);
    else % face = 3, 5, or 7
        y = index_from_coord(x_coord - 1, y_coord + 1, face, N);
    end
end

function w = row_width(index, N) 
    face = face_from_index(index, N);
    y = yfacecoord(index, N);
    if ((face == 1) || (face == 2) || (face == 4) || (face == 6))
        w = 1 + 2 * (y - 1);
    else
        w = N - 2 * (y - 1);
    end
end


% Gives y coordinate of index within each face
function y = yfacecoord(index, N) 
    sites_per_face = ((N + 1) / 2)^2;
    if (index <= sites_per_face)
        y = ceil(sqrt(index));
    elseif (index > 7 * sites_per_face)
        y = (N + 1) / 2 + 1 - ceil(sqrt(2 * (N + 1)^2 + 1 - index));
    else
        % remove bottom triangle
        index = index - sites_per_face;
        row_length = 3 * (N + 1);

        y = floor((index - 1) / row_length) + 1;
    end
end

% Gives face containing given index CHECKED
function f = face_from_index(index, N) 
    sites_per_face = ((N + 1) / 2)^2;
    if (index <= sites_per_face)
        f = 1;
    elseif (index > 7 * sites_per_face)
        f = 8;
    else
        % remove bottom triangle
        index = index - sites_per_face;
        row_length = 3 * (N + 1);

        face_pair = mod(ceil(index / (N + 1)) - 1, 3) + 1;
        row_num = ceil(index / row_length);

        short_width = 1 + 2 * (row_num - 1);
        long_width = (N + 1) - short_width;

        index = mod(index - 1, short_width + long_width) + 1;
        if (index <= short_width)
            f = 2 * face_pair; 
        else
            f = 2 * face_pair + 1;
        end
    end
end

% From face coordinate and face number give index CHECKED
function i = index_from_coord(x, y, face, N)
    sites_per_face = ((N + 1) / 2)^2;
    total_sites = 2 * (N + 1)^2;

    if (face == 1)
        i = (y - 1)^2 + x;
    elseif (face == 2)
        i = sites_per_face + (y - 1) * (3 * (N + 1)) + x;
    elseif (face == 3)
        i = sites_per_face + (y - 1) * (3 * (N + 1)) + (1 + 2 * (y - 1)) + x;
    elseif (face == 4)
        i = sites_per_face + (y - 1) * (3 * (N + 1)) + (N + 1) + x;
    elseif (face == 5)
        i = sites_per_face + (y - 1) * (3 * (N + 1)) + (N + 1) + (1 + 2 * (y - 1)) + x;
    elseif (face == 6)
        i = sites_per_face + (y - 1) * (3 * (N + 1)) + 2 * (N + 1) + x;
    elseif (face == 7)
        i = sites_per_face + (y - 1) * (3 * (N + 1)) + 2 * (N + 1) + (1 + 2 * (y - 1)) + x;
    elseif (face == 8)
        i = total_sites - ((N + 1) / 2 + 1 - y)^2 + x;
    end
end