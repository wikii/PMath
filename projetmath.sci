//Purpose: fchier scilab qui implémente l'algorithme PageRank


//fonction qui renvoie 0 si l'argument est négatif et 1 sinon
function r=r01(t)
    if (t > 0.5) then
        r = 1;
    else r = 0;
    end
endfunction

//Renvoie 0 ou 1 aléatoirement
function r=rand01(t)
    r = r01(rand(t));
endfunction

//fonction qui renvoie un vecteur colone de taille n avec des zéros ou des 1 (aléatoirement)
function r=vectrand(n)
    r = zeros(n)
    for i = 1:n
        r(i)=rand01(n);
    end
endfunction

//fonction qui renvoie une matrice de taille n avec des zéros ou des 1 (aléatoirement)
function r=matricerand(n)
    r = zeros(n,n)
    for i = 1:n
        for j = 1:n
            r(i,j)=rand01(n);
        end
    end
endfunction

//fonction qui renvoie une matrice de la forme de C (aléatoire)
//C est formée de 0 sur sa diagonale et de 0 ou 1 (aléatoire) partout ailleurs
function r=c(n)
    r = matricerand(n)
    for i = 1:n
        r(i,i) = 0;
    end
endfunction

//fonction qui rend la colonne d'une matrice nulle
//Cela permetras de faire des tests
function r=vectnull(M, j)
    n = sqrt(length(M));
    r = M;
    for i = 1:n
        r(i,j)=0;
    end
endfunction

//erreur dans le sujet PAS DE MULTIPLICATION PAR n
function r=Nj2(C, j)
    n = sqrt(length(C));
    r = 0;
    for k = 1:n
        r = r + n*C(k,j)
    end
endfunction

// fonction qui renvoie Nj, le nombre de 1 dans une colonne (=le nombre de liens qui mènent à la page j)
function r=Nj(C, j)
    n = sqrt(length(C));
    r = 0;
    for k = 1:n
        r = r + C(k,j)
    end
endfunction

// fonction qui renvoie Nj, le nombre de 1 dans une colonne (=le nombre de liens qui mènent à la page j) pour un vecteur colone
function r=NjC(C)
    n = length(C);
    r = 0;
    for k = 1:n
        r = r + C(k);
    end
endfunction

// renvoie 0 Si Nj est nul sinon renvoie le coefficient i,j de la matrice C (de départ) divisé par Nj
function r=recNj(C, i, j)
    n = sqrt(length(C));
    f = Nj(C, j);
    if (f==0)
        r=0;
    else r= C(i,j) / Nj(C, j);
    end;
endfunction

// renvoie 0 Si Nj est nul sinon renvoie le coefficient i,j de la matrice C (de départ) divisé par Nj pour un vecteur colonne
function r=recNjC(C, i)
    n = length(C);
    f = NjC(C);
    if (f==0)
        r=0;
    else r= C(i) / NjC(C);
    end;
endfunction

//Renvoie la matrice Q (du sujet) 
function r=q(C)
    n = sqrt(length(C));
    r = zeros(n,n);
    for i = 1:n
        for j=1:n
            r(i,j)=recNj(C,i, j);
        end
    end
endfunction

//Renvoie la matrice Q (du sujet) pour un vecteur colone
function r=qC(C)
    n = length(C);
    r = zeros(n);
    for i = 1:n
        r(i)=recNjC(C,i);
    end
endfunction

//fonction qui renvoie le vecteur colonne e (avec que des 1)
function r=e(n)
    r = ones(n,1);
endfunction

//fonction qui renvoie le vecteur ligne te, la transposée de e (avec que des 1)
function r=te(n)
    r = ones(1,n);
endfunction

//fonction qui renvoie les dj, dj vaut 1 si Nj vaut 0 et il vaut 0 autrement
function r=dj(C, j)
    n = sqrt(length(C));
    f = Nj(C, j);
    if (f==0)
        r=1;
    else r=0;
    end
endfunction

//fonction qui renvoie le vecteur colonne d qui contient les dj de chaque colonne
function r=d(C)
    n = sqrt(length(C));
    r = zeros(1,n);
        for j = 1:n
            r(j) = dj(C, j);
        end
endfunction

//fonction qui renvoie la matrice P: c'est la matrice Q dont les colonnes nulles sont remplaceés par un vecteur colonne 1/n
function r=p(C)
    n = sqrt(length(C));
    r = zeros(n,n);
    Q = q(C);
    D = e(n)*d(C);
    r = Q + (1/n)*D;
endfunction

//fonction qui verifie que la somme des coefficients de chaque colonne de P vaut bien 1
function r=verifp(C)
    n = sqrt(length(C));
    f = 0;
    r = "non";
    P = p(C);
    for i = 1:n
        f = f + sum(P(:,i));
    end
    if (f == n) then
        r = "oui";
    end
endfunction

//fonction suplémentaire (pour faire joli) : fonction qui test si on a bien 1 valeur propre de P
function r= verifvalp(C)
    n = sqrt(length(C));
    r = "faux";
    f = te(n)*p(C);
    if f == te(n) then
        r ="vrai";
    end
endfunction

//fonction qui renvoie la matrice A du sujet
function r=a(C, b)
    n = sqrt(length(C));
    r = zeros(n,n);
    P = p(C);
    E = te(n)*e(n);
    r = b*P + (1-b)*(1/n)*E;
endfunction

//fonction qui renvoie le vecteur propre d'une matrice C associé à la valeur propre l
function r=vectpropre(C, l)
    n = sqrt(length(C));
    r = kernel(C-eye(n,n));
endfunction

//fonction puissance 
// Calcul du vecteur solution r par la méthode des puissances mis dans une matrice 
function [r,it] = methodePuissance(A, ordreGrandeur)
    n = sqrt(length(A));
    x0 = vectrand(n);
    r = x0; //Initialisation, on peut prendre différentes valeurs mais pas le vecteur nul (fonction spécial pour le faire) doit etre formé de zero et de 1
    x = A*r ; //Première itération
    it = 0 ;
    while abs(r-x)>10^(-ordreGrandeur) //Tant que la "convergence" est supérieure à l'ordre de grnadeur souhaité
        r = x; //Stockage de l'itération actuelle
        x = A*x/norm(A*x); //Calcul de la prochaine itération
        it = it + 1;
    end
    r = x //Stockage de la valeur finale
endfunction

//reccuperation du vecteur propre associé à 1
function r=v(A, og)
    n = sqrt(length(A));
    B = methodePuissance(A,og);
    for (i = 1:n)
        r(i) = B(i, 1);
    end
endfunction

function [r,it] = methodePuissanceN(A, ordreGrandeur)
    n = sqrt(length(A));
    x0 = qC(vectrand(n)); //
    r = x0; //Initialisation, on peut prendre différentes valeurs mais pas le vecteur nul (fonction spécial pour le faire) doit etre formé de zero et de 1
    x = A*r ; //Première itération
    it = 0 ;
    while abs(r-x)>10^(-ordreGrandeur) //Tant que la "convergence" est supérieure à l'ordre de grnadeur souhaité
        r = x; //Stockage de l'itération actuelle
        x = A*x/norm(A*x); //Calcul de la prochaine itération
        it = it + 1;
    end
    r = x //Stockage de la valeur finale
endfunction

//reccuperation du vecteur propre associé à 1
function r=vN(A, og)
    n = sqrt(length(A));
    B = methodePuissanceN(A,og);
    for (i = 1:n)
        r(i) = B(i, 1);
    end
endfunction
