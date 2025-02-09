function phi_eval =  eval_basis(index, values, Length)
    % given values = [x1 x2 ... xn]
    % it returns [\phi_i (x1) ... \phi_i (xn)]
    % where i is the index
    % Length is the Lx or Ly used to define \phi. The point is that this
    % works for both x and y by exploiting symmetry
    phi_eval = cos(index * pi * values/Length) * sqrt(2/Length);
    if index == 0
        phi_eval = 1/sqrt(Length) * ones(length(values),1)';
    end
end

