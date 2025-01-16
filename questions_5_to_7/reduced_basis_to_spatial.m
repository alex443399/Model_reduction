function T_profile = reduced_basis_to_spatial(a, basis_in_matrix_form)
    T_profile = zeros(size(squeeze(basis_in_matrix_form(1,:,:))));
    R = length(a);
    for r = 1:R
        T_profile = T_profile + a(r) * squeeze(basis_in_matrix_form(r,:,:));
    end
end