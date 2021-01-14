#sequences a l'etude
seq1 = input("Entrez la premiere sequence: ")
seq2 = input("Entrez la seconde sequence: ")

#longueur de sequences definissent dimensions
hauteur = len(seq1)
largeur = len(seq2)

#les tableaux
tableau = {}
traceback = {}

#initialisation des tableaux
for i in range(hauteur + 1): #+1 pour inclure derniere valeur (dernier index)
    tableau[(i, 0)] = i
    traceback[(i,0)] = "U" #vers le haut

for j in range(largeur + 1):
    tableau[(0, j)] = j
    traceback[(0,j)] = "L" #vers la gauche


for i in range(1, hauteur + 1): # on commence a 1 car la colonne 1 est rempli
    for j in range(1, largeur + 1): #on commence a 1 car la rangee 1 est rempli
        if seq1[i - 1] == seq2[j - 1]:#si les 2 nucleotides sont egales..
             sub = 0
        else:
             sub = 1

        cout_ins_U = tableau[(i - 1, j)] + 1 #case haut
        cout_suppr_L = tableau[(i, j - 1)] + 1 #case gauche
        cout_sub_D = tableau[(i - 1, j - 1)] + sub #case diagonale

        cout_min = min(cout_ins_U, cout_suppr_L, cout_sub_D) #retenir valeur minimale entre les 3 dans position
        tableau[(i, j)] = cout_min #placer valeur min dans (i,j)

        if cout_min == cout_ins_U: #si min = Up ; enregistre "U"
            traceback[(i,j)] = "U"

        elif cout_min == cout_suppr_L: #si min = Left; enregistre "L"
            traceback[(i,j)] = "L"

        elif cout_min == cout_sub_D: #si min = D = Diagonal; enregistre "D"
            traceback[(i,j)] = "D"

#calcul distance
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
