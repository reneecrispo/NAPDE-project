addpath (genpath ('nurbs-1.4.3'))
addpath (genpath ('geopdes'))
addpath (genpath ('geopdes_hierarchical'))

addpath (genpath ('Esercizi_Dede'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% ESERCITAZIONE 2 %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ES1
% nel primo esercizio sono forniti dei knots vectors e chiede la regolarità
% delle B-splines associate. Si risolve semplicemente con la formula che
% lega p e n. Roba teorica.

% ES2 
% per tracciare le basis function definite nell'esercizio precedente
% bisogna usare il file:

% display_univariate_nurbs_basis_function

% semplicemente si impostano i parametri desiderati, crea le basis function
% tramite funzioni della libreria e ne traccia i grafici. 
% DOMANDA: SERVE CAPIRE COME FUNZIONA E COME CREA LE BASI? 

% ES 3 
% applica le formule per trovare le espressioni delle basis functions in
% modo analitico, semplicemente applica una formula. 

% ES 4
% inizia a rappresentare le curve con questo metodo usando la funzione:

% display_nurbs_curve

% come prima, setto i parametri che mi servono (knot vector, punti dello 
% spazio fisico, numero basis function, pesi) e poi fa tutto lei. 
% DOMANDA: in base a cosa sono fissati i punti dello spazio fisico al fine
% di rappresentare le curve richieste? 

% ES 5 
% basis functions con pesi, non cambia nulla 

% ES 6 
% chiede di rappresentare un cerchio -> come es 4 ma con pesi

% ES7 
% chiede di disegnare una spirale in 3D- > come es 6 ma con più punti,
% anche in 3D, per il resto non cambia nulla. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% ESERCITAZIONE 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ES 1
% chiede di rappresentare con nurbs delle superfici. per farlo utilizza il file

% display_nurbs_surface

% come negli altri file, basta settare i parametri per ottenere quello che
% ti serve.

% punto a: 
% note: 
% thick indica l'ampiezza cioè il raggio del settore circolare desiderato (provando a cambairlo)
% nel primo step crea il cerchio, nel secondo espande in direzione radiale
% formando la superficie. 
% punto b è uguale ma in 3D perchè fa la supeficie del cilindro, cambiano i punti e i
% valori di height e thick
% punto c idem, varia solo ttt e h

% DOMANDA: in base a cosa assegno i valori di thick??????

%ES 2
% rappresenta una corona cilindrica (solido) in 3D usando il file
% display_nurbs_solid, a noi non interessa

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% ESERCITAZIONE 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% riguarda i k,h,p refinement, per ora non so se ci serve.
% mi sembrano un po' complicate, se serviranno le guarderò meglio :/


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% ESERCITAZIONE 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% è un po' particolare perchè inizia a trattare con funzioni definite con
% le basis function. 
% l'es 1 secondo me è poco importante, il 2 vale la pena guardarlo perchè
% inizia ad inserire equazioni e matrici, almeno la parte introduttiva. é
% utile perchè scompone la funzione cercata come combinazione lineare delle
% basi con opportuni coefficienti. Per trovare i coefficienti si risolve un
% sistema lineare. (sono solo calcoli su carta, nulla su matlab) 
% Il 3 è applicazione diversa del metodo utilizzato nel 2.


