SIMULAZIONI

results_a: 
    mu = 1e-1;
    u_iniziale = rand();

results_b: 
    mu = 1e1;
    u_iniziale = 0.2 + 0.6 * rand();
    Commenti: tutti i valori erano tra 0 e 1 

results_c: 
    mu = 1e1
    u_iniziale = 0.2 + 0.6 * rand();
    Commenti: newton per sistemi non lineari, valori tra 0.2 e 0.8 