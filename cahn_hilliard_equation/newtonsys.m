function [x_new, it] = newtonsys(x0, nmax, toll, fun, J)
% x0 è un vettore colonna
% J matrice jacobiana

    x_new = x0;
    err0 = norm(fun(x0));
    err = 1;
    
    it = 0;
    x_old = x0;

    while ( it < nmax && ( err > toll && err0 > 1e-8 ) )      
        it = it + 1;        
        x_new = x_old - J(x_old) \ fun(x_old);
        err = norm(fun(x_new)) / err0;
        x_old = x_new;        
    end
end