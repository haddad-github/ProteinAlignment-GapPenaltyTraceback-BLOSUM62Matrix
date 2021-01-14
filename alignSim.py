#sequences a l'etude
seq1 = input("Entrez la premiere sequence: ")
seq2 = input("Entrez la seconde sequence: ")

#longueur de sequences definissent dimensions
hauteur = len(seq1)
largeur = len(seq2)

#les tableaux
tableau = {}
traceback = {}

matrice = {"A":{"A":2, "C":-1, "G":-1, "T":-1}, "C":{"A":-1, "C":2, "G":-1, "T": -1},\
          "G":{"A":-1, "C":-1, "G":2, "T":-1}, "T":{"A":-1, "C":-1, "G":-1, "T":2}}
          #matrice pour chaque lettre -> 4 possibilitees

indel = -2 #cout

#initialisation des tableaux
for i in range(hauteur + 1): #+1 pour inclure derniere valeur (dernier index)
    tableau[(i, 0)] = i * indel #le cout indel pour chaque i dans la colone
    traceback[(i,0)] = "U" #vers le haut

for j in range(largeur + 1):
    tableau[(0, j)] = j * indel #le cout indel pour chaque j dans la rangee
    traceback[(0,j)] = "L" #vers la gauche

#initialisation des tableaux
for i in range(1, hauteur + 1): # on commence a 1 car la colonne 1 est rempli
    for j in range(1, largeur + 1): #on commence a 1 car la rangee 1 est rempli
        sub = matrice[seq1[i-1]][seq2[j-1]] #sub = [position nucleotide sequence 1][" sequence 2]

        cout_ins_U = tableau[(i - 1, j)] + indel #+cout indel
        cout_suppr_L = tableau[(i, j - 1)] + indel #+cout indel
        cout_sub_D = tableau[(i - 1, j - 1)] + sub #case diagonale

        cout_max = max(cout_ins_U, cout_suppr_L, cout_sub_D) #retenir valeur minimale entre les 3 dans position
        tableau[(i, j)] = cout_max #placer valeur min dans (i,j)

        if cout_max == cout_ins_U: #si max = Up ; enregistre "U"
            traceback[(i,j)] = "U"

        elif cout_max == cout_suppr_L: #si max = Left; enregistre "L"
            traceback[(i,j)] = "L"

        elif cout_max == cout_sub_D: #si max = D = Diagonal; enregistre "D"
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
