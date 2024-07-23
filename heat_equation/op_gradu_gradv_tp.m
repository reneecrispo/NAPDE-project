function varargout = op_gradu_gradv_tp (spu, spv, msh)

  gradu_u = reshape(spu.shape_function_gradients, spu.ncomp, [], ...
                    msh.nqn, spu.nsh_max, msh.nel);
  gradu_v = reshape(spv.shape_function_gradients, spv.ncomp, [], ...
                    msh.nqn, spv.nsh_max, msh.nel);

  rows = zeros(msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros(msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros(msh.nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if all(msh.jacdet(:, iel))
      jacdet_weights = reshape(msh.jacdet(:, iel) .* ...
                               msh.quad_weights(:, iel), 1, msh.nqn);

      gradu_u_iel = permute(gradu_u(:, :, :, 1:spu.nsh(iel), iel), [1 2 4 3]);
      gradu_u_iel = reshape(gradu_u_iel, spu.ncomp * size(gradu_u, 2), spu.nsh(iel), msh.nqn);
      gradu_u_iel = permute(gradu_u_iel, [1 3 2]);

      gradu_v_iel = permute(gradu_v(:, :, :, 1:spv.nsh(iel), iel), [1 2 4 3]);
      gradu_v_iel = reshape(gradu_v_iel, spv.ncomp * size(gradu_v, 2), spv.nsh(iel), msh.nqn);
      gradu_v_iel = permute(gradu_v_iel, [1 3 2]);

      for ish = 1:spv.nsh(iel)
        gradv_dot_gradu_iel(:, :, ish) = gradu_u_iel(:, :, ish)' * gradu_v_iel(:, :, ish);
      end

      for idof = 1:spv.nsh(iel)
        rows(ncounter+(1:spu.nsh(iel))) = spv.connectivity(idof, iel);
        cols(ncounter+(1:spu.nsh(iel))) = spu.connectivity(1:spu.nsh(iel), iel);

        aux_val = bsxfun(@times, jacdet_weights, gradv_dot_gradu_iel(:, :, idof));
        values(ncounter+(1:spu.nsh(iel))) = sum(sum(aux_val, 2), 1);
        ncounter = ncounter + spu.nsh(iel);
      end
    else
      warning('geopdes:jacdet_zero_at_quad_node', ...
              'op_gradu_gradv_tp: singular map in element number %d', iel)
    end
  end

  if nargout == 1
    varargout{1} = sparse(rows(1:ncounter), cols(1:ncounter), ...
                          values(1:ncounter), spv.ndof, spu.ndof);
  elseif nargout == 3
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error('op_gradu_gradv_tp: wrong number of output arguments')
  end

end
