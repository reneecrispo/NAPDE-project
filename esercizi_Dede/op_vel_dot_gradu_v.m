% OP_VEL_DOT_GRADU_V: assemble the matrix B = [b(i,j)], b(i,j) = (V * grad u_j, v_i), exploiting the tensor product structure.
%
%   mat = op_vel_dot_gradu_v (spu, spv, msh, Vel);
%   [rows, cols, values] = op_vel_dot_gradu_v (spu, spv, msh, Vel);
%
% INPUT:
%
%   spu:   structure representing the space of trial functions (see sp_bspline_2d/sp_evaluate_col)
%   spv:   structure representing the space of test functions (see sp_bspline_2d/sp_evaluate_col)
%   msh:   structure containing the domain partition and the quadrature rule (see msh_2d/msh_evaluate_col)
%   Vel:   advection field
%
% OUTPUT:
%
%   mat: assemble advection matrix
%   rows:   row indices of the nonzero entries
%   cols:   column indices of the nonzero entries
%   values: values of the nonzero entries
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, Rafael Vazquez
% Copyright (C) 2013, Luca Dede'
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function varargout = op_vel_dot_gradu_v (spu, spv, msh, coeff)

  gradu = reshape (spu.shape_function_gradients, spu.ncomp, [], ...
		   msh.nqn, spu.nsh_max, msh.nel);
  shpv = reshape (spv.shape_functions, spv.ncomp, ...
                  msh.nqn, spv.nsh_max, msh.nel);
  
  ndir = size (gradu, 2);

  rows = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  cols = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);
  values = zeros (msh.nel * spu.nsh_max * spv.nsh_max, 1);

  ncounter = 0;
  for iel = 1:msh.nel
    if (all (msh.jacdet(:, iel)))
      jacdet_weights = reshape (msh.jacdet(:, iel) .* ...
                       msh.quad_weights(:, iel), 1, msh.nqn);
                                      
      gradu_iel = permute (gradu(:, :, :, 1:spu.nsh(iel), iel), [1 2 4 3]);
      gradu_iel = reshape (gradu_iel, spu.ncomp * ndir, spu.nsh(iel), msh.nqn);
      gradu_iel = permute (gradu_iel, [1 3 2]);
      
      coeff_iel = reshape (coeff(:, iel), msh.nqn, ndir );
      coeff_iel = permute (coeff_iel, [2 1]);

      for ish = 1:spu.nsh(iel)
        vel_dot_gradu_iel( :, :, ish ) = coeff_iel .* gradu_iel( :, :, ish );
      end
      
      shpv_iel = reshape (shpv(:, :, 1:spv.nsh(iel), iel), spv.ncomp, msh.nqn, spv.nsh(iel));
      shpv_times_jw = bsxfun (@times, jacdet_weights, shpv_iel);
            
      for idof = 1:spv.nsh(iel)
        rows(ncounter+(1:spu.nsh(iel))) = spv.connectivity(idof, iel);
        cols(ncounter+(1:spu.nsh(iel))) = spu.connectivity(1:spu.nsh(iel), iel);

        aux_val = bsxfun (@times, shpv_times_jw(:,:,idof), vel_dot_gradu_iel);
        values(ncounter+(1:spu.nsh(iel))) = sum (sum (aux_val, 2), 1);
        ncounter = ncounter + spu.nsh(iel);
      end
    else
      warning ('geopdes:jacdet_zero_at_quad_node', 'op_vel_dot_gradu_v: singular map in element number %d', iel)
    end
  end

  if (nargout == 1)
    varargout{1} = sparse (rows(1:ncounter), cols(1:ncounter), ...
                           values(1:ncounter), spv.ndof, spu.ndof);
  elseif (nargout == 3)
    varargout{1} = rows(1:ncounter);
    varargout{2} = cols(1:ncounter);
    varargout{3} = values(1:ncounter);
  else
    error ('op_vel_dot_gradu_v: wrong number of output arguments')
  end

end
