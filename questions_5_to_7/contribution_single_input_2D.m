function contribution_matrix = contribution_single_input_2D( ...
    K, L, X, Y, W, Lx, Ly)
    contributions_x = contributions_single_input_1D( ...
        K, X, W, Lx);
    contributions_y = contributions_single_input_1D( ...
        L, Y, W, Ly);
    contribution_matrix = contributions_y'*contributions_x;
end

function contributions = contributions_single_input_1D( ...
    Truncation_size, center, Width, Length_of_rectangle)
    % Truncation size corresponds to K or L
    % Center to X1, Y1, X2, Y2, width to W, Length_of_rectangle to Lx or Ly

    indeces = 0:Truncation_size;
    indeces_scaled_for_cosine = indeces * pi/Length_of_rectangle;
    cosine_component = cos(indeces_scaled_for_cosine*center);
    sine_component = sin(indeces_scaled_for_cosine*Width/2);
    contributions = cosine_component.*sine_component./indeces;
    contributions = contributions/pi*sqrt(8*Length_of_rectangle);
    
    contributions(1) = Width/sqrt(Length_of_rectangle);
end

