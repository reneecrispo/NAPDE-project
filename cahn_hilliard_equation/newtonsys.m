function [xvect, it] = newtonsys(x0, nmax, toll, fun, J)
% x0 è un vettore colonna
% J matrice jacobiana

    xvect = x0;
    err = toll+1;
    it = 0;
    x_old = x0;

    while it < nmax && err > toll
        if(det(J(x_old))==0)
            error ('Matrice non singolare')
        end
        
        it = it + 1;
        x_new = x_old -J(x_old)\fun(x_old);
        err = abs(max(fun(x_new)));
        xvect = [xvect, x_new];
        x_old = x_new;
    end
    
    fprintf('Numero di Iterazioni: %d\n', it);
    fprintf('Radice calcolata: %-12.8g\n', x_new);
end