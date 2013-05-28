// Définir le vecteur X :
//Sur un exemple concret
//Ou sur :
X = matriceX(5)

// Définition des différents vecteurs et matrices nécessaires
N = sommecolonne(X)
Nb = tailleEspace(X)
Q = matriceQ(X, N)
e = vecteurE(X)
d = vecteurD(N)
P = matriceP(Q, e, d, Nb)
A = matriceA(P, Nb, e, 0.85)

// Définition du vecteur initial x0
x0 = rand(Nb, 1)

// Méthode de la puissance
[r1, it1] = methodePuissance(A, x0, 5)

// Méthode de la puissance avec la première amélioration
[r2, it2] = methodeAmelioree(A, x0, 5)

// Méthode de la puissance avec les 2 améliorations
r3 = methodeAmelioree2(Q, e, x0, 0.85, Nb, 30)
