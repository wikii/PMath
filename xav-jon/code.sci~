// Modélisation du sujet

// Création d'une matrice carrée aléatoire de taille n et de coefficient 0 ou 1
function X = matriceX(n)
    X = rand(n,n);
    X = round(X);
    X = - X .* (eye(n,n)-ones(n,n)) ;
endfunction

// Création d'un vecteur qui somme les valeurs de chaque colonne d'une matrice
function N = sommecolonne(X)
    N = sum(X, 'r');
endfunction

// Taille de l'espace dans lequel on travaille
function Nb = tailleEspace(X)
    Nb = sqrt(length(X));
endfunction

// Création de la matrice des pondérations associé à la matrice X
function Q = matriceQ(X, N)
    n = sqrt(length(X)); //Récupération de la taille de la matrice X 
    for i=1:n
        for j=1:n
            if (N(j)==0) then   //Si la j-ème colonne est nulle
                Q(i,j) = 0; 
            else                //Sinon
                Q(i,j) = X(i,j)/N(j) ;
            end
        end
    end
endfunction

// Création du vecteur e : vecteur de "1" de 
function e = vecteurE(X)
    e = ones(sqrt(length(X)), 1);
endfunction

// Création du vecteur d associé à la matrice X, via le vecteur N
function d = vecteurD(N)
    n = length(N); //Récupération de la taille de N
    for i=1:n
        if (N(i)==0) then //Si la i-ème colonne est nulle
            d(i,1) = 1 ;
        else              //Sinon
            d(i) = 0 ;
        end
    end
endfunction

// Création de la matrice P
function P = matriceP(Q, e, d, Nb)
    P = Q + (1/Nb)*e*d';
endfunction

// Calcul du vecteur z qui vérifie l'exactitude du calcul de la matrice P
function z = verifValeurPropre(P, e)
    z = e'*P-e';
endfunction

// Création de la matrice A
function A = matriceA(P, Nb, e, alpha)
    A = P*alpha + ((1-alpha)/Nb)*e*e';
endfunction

// Fonction qui créé directement la matrice A à partir de la matrice X
function A = test(X)
    N = sommecolonne(X)
    Nb = sommeN(N)
    Q = matriceQ(X, N)
    e = vecteurE(X)
    d = vecteurD(N)
    P = matriceP(Q, e, d, Nb)
    z = verifValeurPropre(P, e)
    A = matriceA(P, Nb, e, 0.85)
endfunction

// Calcul du vecteur solution r par la méthode des puissances
function [r,it] = methodePuissance(A, x0, ordreGrandeur)
    r = x0; //Initialisation
    x = A*r ; //Première itération
    it = 0 ;
    while abs(r-x)>10^(-ordreGrandeur) //Tant que la "convergence" est supérieure à l'ordre de grnadeur souhaité
        r = x; //Stockage de l'itération actuelle
        x = A*x/norm(A*x); //Calcul de la prochaine itération
        it = it + 1;
    end
    r = x //Stockage de la valeur finale
endfunction

// Tache 2

// Question 1

// Calcul du vecteur solution r par la méthode des puissances améliorée
function [r,it] = methodeAmelioree(A, x0, ordreGrandeur)
    r = x0./norm(x0); //Initialisation en normalisant
    x = A*r ; 
    it = 0 ;
    while abs(r-x)>10^(-ordreGrandeur)
        r = x; 
        x = A*x ; //Calcul de la prochaine itération sans renormaliser
        it = it + 1;
    end
    r = x ; //Stockage de la valeur finale
endfunction


// Question 2

// Calcul du vecteur solution r par la méthode des puissances  encore améliorée
function r = methodeAmelioree2(Q, e, x0, alpha, Nb, it)
    r = x0./norm(x0,1);
    tmp = ones(Q) ; //Matrice temporaire
    somme = zeros(Q); //Somme des matrices temporaires
    for i=0:it
        somme = somme + tmp ; //Ajout de l'itération actuelle
        tmp = alpha * Q * tmp; //Calcul de la nouvelle itération
    end
    r = tmp * r + somme * (((1-norm(alpha*Q,1))/Nb).*e) ; //Calcul de la valeur finale
endfunction