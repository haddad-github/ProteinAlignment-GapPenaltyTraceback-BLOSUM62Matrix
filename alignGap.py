import math

#sequences a l'etude
seq1 = input("Entrez la premiere sequence: ")
seq2 = input("Entrez la seconde sequence: ")

#longueur de sequences definissent dimensions
hauteur = len(seq1)
largeur = len(seq2)

#les tableaux
tableau = {}
traceback = {}
tableau_Left = {}
tableau_Up = {}

matrice = {"A":{"A":2, "C":-1, "G":-1, "T":-1}, "C":{"A":-1, "C":2, "G":-1, "T": -1},\
          "G":{"A":-1, "C":-1, "G":2, "T":-1}, "T":{"A":-1, "C":-1, "G":-1, "T":2}}
          #matrice pour chaque lettre -> 4 possibilitees

#couts
cout_Ouverture = -2
cout_Extension = -1

#initialisation des tableaux
for i in range(hauteur + 1): #+1 pour inclure derniere valeur (dernier index)

    w_i = cout_Ouverture + (cout_Extension*(i-1)) #formule w en i

    tableau[(0,0)] = 0 #premiere valeur = 0
    tableau[(i,0)] = w_i #formule w

    traceback[(0,0)] = "X"
    traceback[(i,0)] = "U"

    tableau_Left[(i,0)] = w_i #formule w
    tableau_Up[(i,0)] = -math.inf #les infinis


for j in range(largeur + 1):

    w_j = cout_Ouverture + (cout_Extension * (j - 1)) #formule w en j

    tableau[(0,0)] = 0 #premiere valeur = 0
    tableau[(0,j)] = w_j #formule w

    traceback[(0,j)] = "L"

    tableau_Left[(0,j)] = -math.inf #les infinis
    tableau_Up[(0,j)] = w_j #formule w


for i in range(1, hauteur + 1): # on commence a 1 car la colonne 1 est rempli
    for j in range(1, largeur + 1): #on commence a 1 car la rangee 1 est rempli

        sub = matrice[seq1[i - 1]][seq2[j - 1]] #sub = [position nucleotide sequence 1][" sequence 2]

        #les formules pour les tableaux gauche, sim, diag
        cout_max_L = max(tableau_Left[(i, j - 1)] + cout_Extension, tableau[(i, j - 1)] + cout_Ouverture) #gauche
        cout_max_U = max(tableau_Up[(i - 1, j)] + cout_Extension, tableau[(i - 1, j)] + cout_Ouverture) #haut
        cout_max_D = tableau[(i - 1, j - 1)] + sub #s_ij

        cout_max = max(cout_max_U, cout_max_L, cout_max_D) #max entre les 3

        if cout_max == cout_max_U: #si max = Up ; enregistre "U"
            traceback[(i,j)] = "U"
        elif cout_max == cout_max_L: #si max = Left; enregistre "L"
            traceback[(i,j)] = "L"
        elif cout_max == cout_max_D: #si max = D = Diagonal; enregistre "D"
            traceback[(i,j)] = "D"

        tableau_Left[(i, j)] = cout_max_L #pour la table Gauche, la case prend la valeur du max_Left
        tableau_Up[(i, j)] = cout_max_U #pour la table Haut, la case prend la valeur du max_Up
        tableau[(i, j)] = cout_max #pour le tableau principal, la case prend la valeur max entre tous

        
#distance
distance = tableau[(hauteur, largeur)] #derniere valeur du tableau = score

#alignement
cell = (hauteur, largeur) #derniere cellule

align_seq1 = ""
align_seq2 = ""

#traceback
while cell != (0,0): #arreter a (0,0)

    (i, j) = cell

    if traceback[cell] == "U": #si on va en haut..
        align_seq1 = seq1[i-1] + align_seq1
        align_seq2 = "-" + align_seq2
        cell = (i-1, j)

    elif traceback[cell] == "L": #si on va a gauche..
        align_seq1 = "-" + align_seq1
        align_seq2 = seq2[j-1] + align_seq2
        cell = (i, j-1)

    else: #sinon (on va en diagonale)
        align_seq1 = seq1[i-1] + align_seq1
        align_seq2 = seq2[j-1] + align_seq2
        cell = (i-1, j-1)

#imprimer tableau
tmp = "  - "
for j in range(1, largeur + 1):
    tmp += seq2[j - 1] + " "
print(tmp)

for i in range(hauteur + 1):
    if i == 0:
        tmp = "- "
    else:
        tmp = seq1[i - 1] + " "
    for j in range(largeur + 1):
        tmp += str(tableau[(i, j)]) + " "
    print(tmp)

print("Score:", distance)
print(align_seq1)
print(align_seq2)
