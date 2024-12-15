function A = jump_initial_conditions(K,L,Lx,Ly)
    A = contribution_single_input_2D(K,L,Lx/2,Ly/2, min(Lx,Ly)/2,Lx,Ly);
end

