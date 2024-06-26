%Testing
N = 5;
total = 5 * (N + 1)^2;

fprintf([num2str(index_from_coord(1, 1, 11, N)) '\n'])

%% Stop
for i = 1:total
    %if (has_upper(i, N) == false)
    
        fprintf([num2str(i) ' ' num2str(c5_index(i, N)) '\n'])
    %end
end

%% Functions

% Gives width of row within face CHECKED
function w = row_width(index, N) 
    face = face_from_index(index, N);
    y = yfacecoord(index, N);
    if ((face == 1) || (face == 2) || (face == 3) || (face == 4) || (face == 5) || (face == 7)...
            || (face == 9) || (face == 11) || (face == 13) || (face == 15))
        w = 1 + 2 * (y - 1);
    else
        w = N - 2 * (y - 1);
    end
end

% Gives y coordinate of index within each face CHECKED
function y = yfacecoord(index, N) 
    sites_per_face = ((N + 1) / 2)^2;
    face = face_from_index(index, N);
    if (index <= 5 * sites_per_face) % bottom row
        index = mod(index - 1, sites_per_face) + 1;
        y = ceil(sqrt(index));
    elseif (index > 15 * sites_per_face)
        total_sites = 5 * (N + 1)^2 - (20 - face) * sites_per_face;
        y = (N + 1) / 2 + 1 - ceil(sqrt(total_sites + 1 - index));
    else
        % remove bottom row of triangle
        index = index - 5 * sites_per_face;
        row_length = 5 * (N + 1);

        y = floor((index - 1) / row_length) + 1;
    end
end

% Gives x coordinate of index within each face CHECKED
function x = xfacecoord(index, N) 
    face = face_from_index(index, N);
    y_coord = yfacecoord(index, N);
    sites_per_face = ((N + 1) / 2)^2;
    row_w = row_width(index, N);

    if (face <= 5) % bottom row
        index = mod(index - 1, sites_per_face) + 1;
        x = index - (y_coord - 1)^2; 
    elseif (face > 15) %top row
        total_sites = 5 * (N + 1)^2 - (20 - face) * sites_per_face;
        x = row_w - (total_sites + 1 - index - ((N + 1) / 2 - y_coord)^2) + 1;
    elseif (face == 6)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1);
    elseif (face == 8)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - (N + 1);
    elseif (face == 10)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - 2 * (N + 1);
    elseif (face == 12)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - 3 * (N + 1);
    elseif (face == 14)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - 4 * (N + 1);
    elseif (face == 7)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - (N + 1 - row_w);
    elseif (face == 9)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - (N + 1 - row_w) - (N + 1);
    elseif (face == 11)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - (N + 1 - row_w) - 2 * (N + 1);
    elseif (face == 13)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - (N + 1 - row_w) - 3 * (N + 1);
    elseif (face == 15)
        x = index - 5 * sites_per_face - (y_coord - 1) * 5 * (N + 1) - (N + 1 - row_w) - 4 * (N + 1);
    end
end

% Gives face containing given index CHECKED
function f = face_from_index(index, N) 
    sites_per_face = ((N + 1) / 2)^2;
    if (index <= 5 * sites_per_face) % bottom row of triangles
        f = ceil(index / sites_per_face);
    elseif (index > 15 * sites_per_face) % top row of triangles
        f = ceil(index / sites_per_face);
    else % in middle row of triangles
        % remove bottom row of triangle
        index = index - 5 * sites_per_face;
        row_length = 5 * (N + 1);

        face_pair = mod(ceil(index / (N + 1)) - 1, 5) + 1;
        row_num = ceil(index / row_length);

        short_width = 1 + 2 * (row_num - 1);
        long_width = (N + 1) - short_width;

        index = mod(index - 1, short_width + long_width) + 1;
        if (index <= long_width)
            f = 2 * face_pair + 4; 
        else
            f = 2 * face_pair + 5;
        end
    end
end

% From face coordinate and face number give index CHECKED
function i = index_from_coord(x, y, face, N)
    sites_per_face = ((N + 1) / 2)^2;

    if (face <= 5)
        i = sites_per_face * (face - 1) + (y - 1)^2 + x;
    elseif (face == 6)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + x;
    elseif (face == 7)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + (N + 1 - (1 + 2 * (y - 1))) + x;
    elseif (face == 8)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + (N + 1) + x;
    elseif (face == 9)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + (N + 1) + (N + 1 - (1 + 2 * (y - 1))) + x;
    elseif (face == 10)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + 2 * (N + 1) + x;
    elseif (face == 11)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + 2 * (N + 1) + (N + 1 - (1 + 2 * (y - 1))) + x;
    elseif (face == 12)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + 3 * (N + 1) + x;
    elseif (face == 13)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + 3 * (N + 1) + (N + 1 - (1 + 2 * (y - 1))) + x;
    elseif (face == 14)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + 4 * (N + 1) + x;
    elseif (face == 15)
        i = 5 * sites_per_face + (y - 1) * (5 * (N + 1)) + 4 * (N + 1) + (N + 1 - (1 + 2 * (y - 1))) + x;
    elseif (face > 15)
        total_sites = 5 * (N + 1)^2 - (20 - face) * sites_per_face;
        i = total_sites - ((N + 1) / 2 + 1 - y)^2 + x;
    end
end


% Functions to find indices of H matrix 
function y = upper(index, N) % CHECKED
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face <= 5 && y_coord < (N + 1) / 2) % in faces 1-5 but not top row
        y = index_from_coord(x_coord + 1, y_coord + 1, face, N);
    elseif (face == 1 && y_coord == (N + 1) / 2) % in face 1 top row
        y = index_from_coord(x_coord, 1, 6, N);
    elseif (face == 2 && y_coord == (N + 1) / 2) % in face 2 top row
        y = index_from_coord(x_coord, 1, 8, N);
    elseif (face == 3 && y_coord == (N + 1) / 2) % in face 3 top row
        y = index_from_coord(x_coord, 1, 10, N);
    elseif (face == 4 && y_coord == (N + 1) / 2) % in face 4 top row
        y = index_from_coord(x_coord, 1, 12, N);
    elseif (face == 5 && y_coord == (N + 1) / 2) % in face 5 top row
        y = index_from_coord(x_coord, 1, 14, N);
    elseif (face > 15) % faces 15-20
        y = index_from_coord(x_coord - 1, y_coord + 1, face, N);
    elseif (face == 7 && y_coord == (N + 1) / 2) % face 7 top row
        y = index_from_coord(x_coord, 1, 16, N);
    elseif (face == 9 && y_coord == (N + 1) / 2) % face 9 top row
        y = index_from_coord(x_coord, 1, 17, N);
    elseif (face == 11 && y_coord == (N + 1) / 2) % face 11 top row
        y = index_from_coord(x_coord, 1, 18, N);
    elseif (face == 13 && y_coord == (N + 1) / 2) % face 13 top row
        y = index_from_coord(x_coord, 1, 19, N);
    elseif (face == 15 && y_coord == (N + 1) / 2) % face 15 top row
        y = index_from_coord(x_coord, 1, 20, N);
    elseif (mod(face, 2) == 0) % face = 6, 8, 10, 12, 14 
        y = index_from_coord(x_coord - 1, y_coord + 1, face, N);
    else % face = 7, 9, 11, 13, 15 not rop row
        y = index_from_coord(x_coord + 1, y_coord + 1, face, N);
    end
end

function y = lower(index, N) % CHECKED
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face <= 5)
        y = index_from_coord(x_coord - 1, y_coord - 1, face, N);
    elseif (face == 6 && y_coord == 1) % face 6 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 1, N);
    elseif (face == 8 && y_coord == 1) % face 8 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 2, N);
    elseif (face == 10 && y_coord == 1) % face 10 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 3, N);
    elseif (face == 12 && y_coord == 1) % face 12 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 4, N);
    elseif (face == 14 && y_coord == 1) % face 14 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 5, N);
    elseif (face == 16 && y_coord == 1) % face 16 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 7, N);
    elseif (face == 17 && y_coord == 1) % face 17 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 9, N);
    elseif (face == 18 && y_coord == 1) % face 18 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 11, N);
    elseif (face == 19 && y_coord == 1) % face 19 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 13, N);
    elseif (face == 20 && y_coord == 1) % face 20 bottom row
        y = index_from_coord(x_coord, (N + 1) / 2, 15, N);
    elseif (face > 15)
        y = index_from_coord(x_coord + 1, y_coord - 1, face, N);
    else
        y = index - 5 * (N + 1) + 1;
    end
end

function y = right(index, N) % CHECKED
    face = face_from_index(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);
    row_w = row_width(index, N);

    if (face == 1 && x_coord == row_w) % face 1 right side
        y = index_from_coord(1, y_coord, 2, N);
    elseif (face == 2 && x_coord == row_w) % face 2 right side
        y = index_from_coord(1, y_coord, 3, N);
    elseif (face == 3 && x_coord == row_w) % face 3 right side
        y = index_from_coord(1, y_coord, 4, N);
    elseif (face == 4 && x_coord == row_w) % face 4 right side
        y = index_from_coord(1, y_coord, 5, N);
    elseif (face == 5 && x_coord == row_w) % face 5 right side
        y = index_from_coord(1, y_coord, 1, N);
    elseif (face == 16 && x_coord == row_w) % face 16 right side
        y = index_from_coord(1, y_coord, 17, N);
    elseif (face == 17 && x_coord == row_w) % face 17 right side
        y = index_from_coord(1, y_coord, 18, N);
    elseif (face == 18 && x_coord == row_w) % face 18 right side
        y = index_from_coord(1, y_coord, 19, N);
    elseif (face == 19 && x_coord == row_w) % face 19 right side
        y = index_from_coord(1, y_coord, 20, N);
    elseif (face == 20 && x_coord == row_w) % face 20 right side
        y = index_from_coord(1, y_coord, 16, N);
    elseif (face == 15 && x_coord == row_w) % face 7 right side
        y = index_from_coord(1, y_coord, 6, N);
    else
        y = index + 1; 
    end
end

function y = left(index, N) % CHECKED
    face = face_from_index(index, N);
    row_w = row_width(index, N);
    x_coord = xfacecoord(index, N);
    y_coord = yfacecoord(index, N);

    if (face == 1 && x_coord == 1) % face 1 left side
        y = index_from_coord(row_w, y_coord, 5, N);
    elseif (face == 2 && x_coord == 1) % face 2 left side
        y = index_from_coord(row_w, y_coord, 1, N);
    elseif (face == 3 && x_coord == 1) % face 3 left side
        y = index_from_coord(row_w, y_coord, 2, N);
    elseif (face == 4 && x_coord == 1) % face 4 left side
        y = index_from_coord(row_w, y_coord, 3, N);
    elseif (face == 5 && x_coord == 1) % face 5 left side
        y = index_from_coord(row_w, y_coord, 4, N);
    elseif (face == 16 && x_coord == 1) % face 16 left side
        y = index_from_coord(row_w, y_coord, 20, N);
    elseif (face == 17 && x_coord == 1) % face 17 left side
        y = index_from_coord(row_w, y_coord, 16, N);
    elseif (face == 18 && x_coord == 1) % face 18 left side
        y = index_from_coord(row_w, y_coord, 17, N);
    elseif (face == 19 && x_coord == 1) % face 19 left side
        y = index_from_coord(row_w, y_coord, 18, N);
    elseif (face == 20 && x_coord == 1) % face 20 left side
        y = index_from_coord(row_w, y_coord, 19, N);
    elseif (face == 6 && x_coord == 1)
        y = index_from_coord(N + 1 - row_w, y_coord, 15, N);
    else
        y = index - 1;
    end  
end

% Indicate whether nearest vertical neighbor is above or below CHECKED
function answer = has_upper(index, N)
    y = yfacecoord(index, N);
    face = face_from_index(index, N);
    sites_per_face = ((N + 1) / 2)^2;

    if (face <= 5) % in lower row 
        index = mod(index - 1, sites_per_face) + 1; 
        if (mod(y, 2) - mod(index, 2) == 0)
            answer = true;
        else 
            answer = false;
        end
    elseif (face > 15) % in upper row
        index = index - (face - 1) * sites_per_face;
        if (mod(y, 2) - mod(index, 2) == 0)
            answer = false;
        else 
            answer = true;
        end
    else % in middle row
        index = index - 5 * sites_per_face;
        if (mod(index, 2) == 0)
            answer = true;
        else
            answer = false;
        end
    end
end

% Function to calculate mean level spacing ratio