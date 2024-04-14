%% Plot ten tetrahedron eigenvectors

% Plot tetrahedron eigenvectors
N = 41;

X = zeros((N + 1)^2, 1);
Y = zeros((N + 1)^2, 1);

a = 1 / N;

% Compute (x,y) coordinates
for i = 1: 2 * (N + 1)^2
    face = face_from_index(i, N);
    x_coord = xfacecoord(i, N);
    y_coord = yfacecoord(i, N);

    x_spacing = 1 / (N + 1);
    y_spacing = sqrt(3) / (N + 1);
    
    if (face == 1)
        X(i, 1) = - (y_coord - 1) * x_spacing + (x_coord - 1) * x_spacing;
        Y(i, 1) = y_spacing / 2 + (y_coord - 1) * y_spacing;
    elseif (face == 2)
        X(i, 1) = - (3/ 2) - (y_coord - 1) * x_spacing + (x_coord - 1) * x_spacing;
        Y(i, 1) = + sqrt(3) / 2 + y_spacing / 2 + (y_coord - 1) * y_spacing;
    elseif (face == 3)
        X(i, 1) = - 1 - (1 / 2) + x_spacing + (y_coord - 1) * x_spacing + (x_coord - 1) * x_spacing;
        Y(i, 1) = sqrt(3) / 2 + y_spacing / 2 + (y_coord - 1) * y_spacing;
    elseif (face == 4)
        X(i, 1) = - (1/ 2) - (y_coord - 1) * x_spacing + (x_coord - 1) * x_spacing;
        Y(i, 1) = + sqrt(3) / 2 + y_spacing / 2 + (y_coord - 1) * y_spacing;
    elseif (face == 5)
        X(i, 1) = - (1 / 2) + x_spacing + (y_coord - 1) * x_spacing + (x_coord - 1) * x_spacing;
        Y(i, 1) = sqrt(3) / 2 + y_spacing / 2 + (y_coord - 1) * y_spacing;
    elseif (face == 6)
        X(i, 1) = 1 / 2 - (y_coord - 1) * x_spacing + (x_coord - 1) * x_spacing;
        Y(i, 1) = + sqrt(3) / 2 + y_spacing / 2 + (y_coord - 1) * y_spacing;
    elseif (face == 7)
        X(i, 1) = 1 - (1 / 2) + x_spacing + (y_coord - 1) * x_spacing + (x_coord - 1) * x_spacing;
        Y(i, 1) = sqrt(3) / 2 + y_spacing / 2 + (y_coord - 1) * y_spacing;
    elseif (face == 8)
        X(i, 1) = - 1 / 2 - 1 / 2 + x_spacing + (y_coord - 1) * x_spacing + (x_coord - 1) * x_spacing;
        Y(i, 1) = sqrt(3) + y_spacing / 2 + (y_coord - 1) * y_spacing;
    end

end

%% Plot the eigenvector
folderName = 'octa_sigma_d_evecs';
mkdir(folderName);

for k = 0:1
    start_index = k * 10 + 1;
    end_index = start_index + 9;
    
    
    for i=start_index:end_index
        cn = 100;          % number of colors
        cm = colormap(parula(cn)); 
        min_of_evec = min(eigenvectors(:, i));
        max_of_evec = max(eigenvectors(:, i));
        color_indices = ceil((max_of_evec - eigenvectors(:, i)) / (max_of_evec - min_of_evec) * cn);
    
        % Remove zeros from color indices
        for j = 1:(2 * (N + 1)^2)
            if (color_indices(j) == 0)
                color_indices(j) = 1;
            end
        end
    
        if (i == start_index)
            figure(1)
            tiledlayout(2, ceil((end_index - start_index + 1)/2), 'TileSpacing', 'tight','Padding','Tight')
            nexttile
            scatter3(X, Y, eigenvectors(:, i), [], cm(color_indices, :), '.')
            xlabel('x');
            ylabel('y');
            view(0,90);
            title(['wavefunction index = ' num2str(i)])
            daspect([1 1 0.2])
        else
            nexttile
            scatter3(X, Y, eigenvectors(:, i), [], cm(color_indices, :), '.')
            xlabel('x');
            ylabel('y');
            view(0,90);
            title(['wavefunction index = ' num2str(i)])
            daspect([1 1 0.2])
        end
    end
    set(figure(1),'position',[0,0,1600,600])
    saveas(gcf, [folderName '/tetra_' num2str(k * 10 + 1) '-' num2str(k * 10 + 10) '.jpeg']);
end

%% Functions
% Gives width of row within face CHECKED
function w = row_width(index, N) 
    face = face_from_index(index, N);
    y = yfacecoord(index, N);
    if ((face == 1) || (face == 2) || (face == 4) || (face == 6))
        w = 1 + 2 * (y - 1);
    else
        w = N - 2 * (y - 1);
    end
end

% Gives y coordinate of index within each face CHECKED
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

% Indicate whether nearest vertical neighbor is above or below CHECKED
function answer = has_upper(index, N)
    y = yfacecoord(index, N);
    face = face_from_index(index, N);
    sites_per_face = ((N + 1) / 2)^2;

    if (face == 1)
        if (mod(y, 2) - mod(index, 2) == 0)
            answer = true;
        else 
            answer = false;
        end
    elseif (mod(sites_per_face, 2) > 0)
        if (face == 8)
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
    else
        if (face == 8)
            if (mod(y, 2) - mod(index, 2) == 0)
                answer = false;
            else 
                answer = true;
            end
        elseif (mod(index, 2) == 0)
            answer = false;
        else
            answer = true;
        end
    end

end


