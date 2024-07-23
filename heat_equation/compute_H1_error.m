function [H1_error, L2_error] = compute_H1_error(space, mesh, u, u_ex, grad_u_ex)
    % Calcola il gradiente della soluzione approssimata u
    [grad_u_x, grad_u_y] = geopdes_grad(u, 'X', mesh, space);

    % Calcola l'errore nel gradiente in norma L2
    grad_error_squared = sum(sum((grad_u_x - grad_u_ex(:,:,1)).^2 + (grad_u_y - grad_u_ex(:,:,2)).^2));

    % Calcola l'errore in norma L2
    L2_error_squared = sum(sum((u - u_ex).^2));

    % Calcolo dell'errore in norma H1
    H1_error = sqrt(L2_error_squared + grad_error_squared);
    L2_error = L2_error_squared;
end