function varargout = op_u_v_tp_cahn_hilliard_non_lin(space1, space2, msh, non_linear_term,u_old)

for idim = 1:msh.ndim
    size1 = size (space1.sp_univ(idim).connectivity);
    size2 = size (space2.sp_univ(idim).connectivity);
    if (size1(2) ~= size2(2) || size1(2) ~= msh.nel_dir(idim))
      error ('One of the discrete spaces is not associated to the mesh')
    end
end
  
A = spalloc (space2.ndof, space1.ndof, 3*space1.ndof);

 %spazio in cui è immerso il problema
  
for iel = 1:msh.nel_dir(1)
   msh_col = msh_evaluate_col(msh, iel);
   sp_col = sp_evaluate_col(space1, msh_col);%per noi trail=test     
     
   coefs = [];
   for l = 1 : msh.nel_dir(2)
     ieloc = iel + msh.nel_dir(2) * ( l - 1 );
     fun_indexes=sp_get_basis_functions(space1,msh,ieloc); %indici delle funzioni 
     %di base che agiscono su quella cella
     %associate alle funzioni di base all'elem iel
     % Estrai i valori di u_old associati agli indici fun_indexes
     u_old_new = u_old(fun_indexes);
     %valutazione delle funzioni di base sui nodi di quadratura     
     ris=sp_col.shape_functions(:, :, iel)'.*u_old_new; %9x9
          
     coefs = [ coefs; sum(ris,1) ];
   end
  coefs = coefs';

  %size_mat=size(shape_fun_for_cells);
   %a questa matrice dovrei applicare il termine non lineare
   for i = 1:msh_col.nqn
    for j = 1:msh_col.nel
         coefs(i, j)=non_linear_term(coefs(i,j));            
    end
   end
  
   
  A = A + op_u_v (sp_col, sp_col, msh_col, coefs);
   
end
 
  if (nargout == 1)
    varargout{1} = A;
  elseif (nargout == 3)
    [rows, cols, vals] = find (A);
    varargout{1} = rows;
    varargout{2} = cols;
    varargout{3} = vals;
  end

end
