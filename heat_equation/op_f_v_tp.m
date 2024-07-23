function rhs_F = op_f_v_tp(space, msh, f)
    % Inizializza il vettore di carico rhs_F
    rhs_F = zeros(space.ndof, 1);
    
    % Itera su tutti gli elementi della mesh
    for iel = 1:msh.nel
        % Ottieni i nodi di quadratura e i pesi
        qn = msh.quad_nodes(:, :, iel);
        qw = msh.quad_weights(:, iel);
        
        % Valuta la funzione sorgente sui nodi di quadratura
        f_val = f(qn(1, :), qn(2, :), 0); % Supponendo che la funzione sorgente non dipenda dal tempo
        
        % Valuta le funzioni di base delle NURBS nei nodi di quadratura
        shp_nurbs = reshape(sp_eval_u(sp_bspline_2d(nrb2bspline(nurbs)), msh, qn), sp.ncomp, msh.nqn, sp.nsh_max, msh.nel);
        
        % Calcola la contribuzione della funzione sorgente all'integrale per l'elemento corrente
        f_contrib = f_val .* repmat(qw', size(f_val, 1), 1);
        
        % Assembla il contributo nel vettore di carico rhs_F
        for idof = 1:space.nsh(iel)
            for j = 1:space.ncomp
                % Ottieni gli indici globali delle funzioni base di test
                row_idx = space.connectivity(idof, iel) + (j - 1) * space.ndof;
                
                % Aggiungi il contributo del nodo di quadratura all'elemento del vettore di carico
                rhs_F(row_idx) = rhs_F(row_idx) + sum(f_contrib(j, :) .* shp_nurbs(j, :, idof, iel));
            end
        end
    end
end
