function [x_vect,it] = newton(x0,nmax,toll,fun,dfun,mol)
    % [x_vect,it] = newton(x0,nmax,toll,fun,dfun,mol)
    
    % Metodo di Newton per la ricerca degli zeri della funzione fun.
    % Test d'arresto basato sull'incremento tra due iterazioni successive.
    
    % Parametri di ingresso:
    
    % x0 Punto di partenza
    % nmax Numero massimo di iterazioni
    % toll Tolleranza sul test d'arresto
    % fun dfun Anonymous function contenenti la funzione e la sua derivata
    % mol Molteplicita' dello zero che si vuole trovare
    
    % Parametri di uscita:
    
    % x_vect Vett. contenente tutte le iterate calcolate
    % (l'ultima componente e' la soluzione)
    % it Iterazioni effettuate

    x_vect = x0;
    err = toll+1;
    it = 0;
    x_old = x0;

    while it < nmax && err > toll
        fx = fun(x_old);
        dfx = dfun(x_old);
        
        if dfx == 0
            % error fa terminare l'esecuzione della function
            error(' Arresto per azzeramento di dfun');
        end

        % nuova iterazione di newton
        x_new = x_old-mol*fx/dfx;
        err = abs(x_new-x_old);

        % salvo le quantita` di interesse negli output
        x_vect = [x_vect; x_new];
        it = it+1;
        x_old = x_new;
    end
    
    fprintf('Numero di Iterazioni: %d\n', it);
    fprintf('Radice calcolata: %-12.8g\n', x_new);
return