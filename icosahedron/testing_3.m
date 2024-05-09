% testing
N = 5;

%x_indexes = [1, 1, 2, 3, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 1];
%y_indexes = [1, 2, 2, 2, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 2, 2, 2, 3, 1, 2, 2, 2, 3, 3, 3, 3, 3, 1, 1, 1, 1, 1, 2, 2, 2, 3];
%faces = [1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4];

x_indexes = [1, 1, 2, 3, 1, 2, 3, 1, 1, 1, 2, 3, 1, 2, 3, 1];
y_indexes = [1, 2, 2, 2, 1, 1, 1, 2, 1, 2, 2, 2, 1, 1, 1, 2];
faces = [1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4];

for i = 1:((N + 1)^2)
    %if (has_upper(i, N) == false)
        fprintf([num2str(i) ' ' mat2str(on_left(i, N)) '\n'])
    %end
end

% Indicate whether index is on left edge
function answer = on_left(index, N)
    x = xfacecoord(index, N);
    face = face_from_index(index, N);

    if (face == 3)
        answer = false;
    elseif (x == 1)
        answer = true;
    else
        answer = false;
    end
end

% Indicate whether index is on right edge
function answer = on_right(index, N)
    x = xfacecoord(index, N);
    row_w = row_width(index, N);
    face = face_from_index(index, N);

    if (face == 2)
        answer = false;
    elseif (x == row_w)
        answer = true;
    else
        answer = false;
    end
end

function i = sigma_d_legit_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);
    if (face == 1) 
        new_face = 1;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 3;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 2;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 4;
        new_x = row_w + 1 - x_coord;
        new_y = y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Given index on discretized octahedron returns new index after inv CHECKED
function i = inv_index(index, N)
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    if (face == 1) 
        new_face = 4;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 2) 
        new_face = 3;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 3) 
        new_face = 2;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    elseif (face == 4) 
        new_face = 1;
        new_x = x_coord;
        new_y = (N + 1) / 2 + 1 - y_coord;
        i = index_from_coord(new_x, new_y, new_face, N);
    end
end

% Indicate whether nearest vertical neighbor is above or below
function answer = has_upper(index, N)
    y = yfacecoord(index, N);
    face = face_from_index(index, N);
    sites_per_face = ((N + 1) / 2)^2;

    if (face == 1 || face == 3) % in lower row 
        index = mod(index - 1, sites_per_face) + 1; 
        if (mod(y, 2) - mod(index, 2) == 0)
            answer = true;
        else 
            answer = false;
        end
    else
        index = index - (face - 1) * sites_per_face;
        if (mod(y, 2) - mod(index, 2) == 0)
            answer = false;
        else 
            answer = true;
        end
    end
end

% CHECKED
function y = left(index, N) 
    face = face_from_index(index, N);
    row_w = row_width(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face == 1 && x_coord == 1) % face 1 left side
        y = index_from_coord(row_w, y_coord, 1, N);
    elseif (face == 2 && x_coord == 1) % face 2 left side
        y = index_from_coord(N + 1 - row_w, y_coord, 3, N);
    elseif (face == 3 && x_coord == 1) % face 3 left side
        y = index_from_coord(N + 1 - row_w, y_coord, 2, N);
    elseif (face == 4 && x_coord == 1) % face 4 left side
        y = index_from_coord(row_w, y_coord, 4, N);
    else
        y = index - 1;
    end  
end

% CHECKED
function y = right(index, N) 
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);

    if (face == 1 && x_coord == row_w) % face 1 right side
        y = index_from_coord(1, y_coord, 1, N);
    elseif (face == 2 && x_coord == row_w) % face 2 right side
        y = index_from_coord(1, y_coord, 3, N);
    elseif (face == 3 && x_coord == row_w) % face 3 right side
        y = index_from_coord(1, y_coord, 2, N);
    elseif (face == 4 && x_coord == row_w) % face 4 right side
        y = index_from_coord(1, y_coord, 4, N);
    else
        y = index + 1; 
    end
end

% Functions to find indices of H matrix CHECKED
function y = upper(index, N) 
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face == 1 && y_coord < (N + 1) / 2) % in face 1 but not top row
        y = index_from_coord(x_coord + 1, y_coord + 1, 1, N);
    elseif (face == 1 && y_coord == (N + 1) / 2) % in face 1 top row
        y = index_from_coord(x_coord, 1, 2, N);
    elseif (face == 3 && y_coord < (N + 1) / 2) % in face 3 but not top row
        y = index_from_coord(x_coord + 1, y_coord + 1, 3, N);
    elseif (face == 3 && y_coord == (N + 1) / 2) % in face 3 top row
        y = index_from_coord(x_coord, 1, 4, N);
    else % faces 2, 4
        y = index_from_coord(x_coord - 1, y_coord + 1, face, N);
    end
end

% From face coordinate and face number give index CHECKED
function i = index_from_coord(x, y, face, N)
    sites_per_face = ((N + 1) / 2)^2;

    if (face == 1 || face == 3)
        i = sites_per_face * (face - 1) + (y - 1)^2 + x;
    else
        total_sites = face * sites_per_face;
        i = total_sites - ((N + 1) / 2 + 1 - y)^2 + x;
    end
end


% Gives x coordinate of index within each face CHECKED
function x = xfacecoord(index, N) 
    face = face_from_index(index, N);
    y_coord = yfacecoord(index, N);
    sites_per_face = ((N + 1) / 2)^2;
    row_w = row_width(index, N);

    if (face == 1 || face == 3) % bottom row
        index = mod(index - 1, sites_per_face) + 1;
        x = index - (y_coord - 1)^2; 
    elseif (face == 2)
        total_sites = 2 * sites_per_face;
        x = row_w - (total_sites + 1 - index - ((N + 1) / 2 - y_coord)^2) + 1;
    elseif (face == 4)
        total_sites = 4 * sites_per_face;
        x = row_w - (total_sites + 1 - index - ((N + 1) / 2 - y_coord)^2) + 1;
    end
end

% Gives y coordinate of index within each face CHECKED
function y = yfacecoord(index, N) 
    sites_per_face = ((N + 1) / 2)^2;
    face = face_from_index(index, N);
    if (face == 1 || face == 3)
        index = mod(index - 1, sites_per_face) + 1;
        y = ceil(sqrt(index));
    else
        index = mod(index - 1, sites_per_face) + 1;
        y = (N + 1) / 2 + 1 - ceil(sqrt(sites_per_face + 1 - index));
    end
end

% Gives face containing given index CHECKED
function f = face_from_index(index, N) 
    sites_per_face = ((N + 1) / 2)^2;
    if (index <= 1 * sites_per_face)
        f = 1;
    elseif (index <= 2 * sites_per_face)
        f = 2;
    elseif (index <= 3 * sites_per_face)
        f = 3;
    else
        f = 4;
    end
end

% Gives width of row within face CHECKED
function w = row_width(index, N) 
    face = face_from_index(index, N);
    y = yfacecoord(index, N);
    if ((face == 1) || (face == 3))
        w = 1 + 2 * (y - 1);
    else
        w = N - 2 * (y - 1);
    end
end