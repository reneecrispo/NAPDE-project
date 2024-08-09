function [xvect, it] = newtonsys_cahn_hilliard(x0, nmax, toll, fun, J)
    % x0 è un vettore colonna
    % J matrice jacobiana (funzione che restituisce una matrice sparsa)

    xvect = x0;
    err = toll+1;
    it = 0;
    x_old = x0;

    while it < nmax && err > toll
        % Calcola la Jacobiana in x_old
        Jx = J(x_old);
        
        % Fattorizzazione LU con pivoting
        [L, U, P] = lu(Jx);

        % Verifica se la matrice è singolare controllando gli elementi diagonali di U
        if any(abs(diag(U)) < eps)
            error('Matrice Jacobiana singolare');
        end
        
        it = it + 1;
        
        % Risolvi il sistema J(x_old) * delta_x = -fun(x_old)
        delta_x = U \ (L \ (P * -fun(x_old)));
        x_new = x_old + delta_x;
        
        err = norm(delta_x, inf); % usa la norma infinito del vettore delta_x per l'errore
        xvect = [xvect, x_new];
        x_old = x_new;
    end
    
    % fprintf('Numero di Iterazioni: %d\n', it);
    % fprintf('Radice calcolata: %-12.8g\n', x_new);
end
