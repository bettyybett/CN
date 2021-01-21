import numpy as np
import math

def SubsDesc(A, b, tol):
    '''
    SubsDesc: rezolva sisteme superior triunghiulare
    Input: a - matricea initiala a sistemului
           b - vectorul termenilor liberi
           tol - valoare cu care vom compara numerele diferite de 0
    output: x -solutia sistemului Ax = b
    '''

    lin = np.shape(A)[0]
    col = np.shape(A)[1]
    if (lin != col):
        print('Matricea nu este patratica, dati alta matrice')
        return None

    # for i in range(lin):
    #     for j in range(i):
    #         if abs(A[i][j]) > tol:
    #             print('Matricea nu este superior triunghiulara, dati o alta matrice')
    #             return None

    for k in range(lin):
        if (abs(A[k][k]) < tol):
            print('Sistemul nu admite solutie unica')
            return None
    # calculam x ul nostru din Ax=b
    x = np.zeros((lin, 1))
    x[lin - 1] = b[lin - 1] / A[lin - 1][lin - 1] #calculam ultima valoare a lui x
    k = lin - 2 ##valoare de la care se calculeaza x ul in continuare
    while (k > - 1):
        sum = 0
        for j in range(k + 1, lin):
            sum += A[k][j] * x[j]
        x[k] = (b[k] - sum) / A[k][k]
        k -= 1
    return x


def GaussPivPart(A, b, tol):
    '''
    GaussFP rezolva sisteme de ecuatii liniare folosind metoda Gauss fara pivotare
    input: A - matricea asociata sistemului(Matrice patratica)
           b - vectorul termenilor liberi
           tol - valoare numerica folosita pentru compararea numerelor apropiate de 0
    output: x - solutia sistemului(Ax = b)
    '''
    # verificam daca matricea este patratica
    lin = np.shape(A)[0]
    col = np.shape(A)[1]
    if (lin != col):
        print('Matricea nu este patratica, dati alta matrice')
        return None, None

    Aext = np.concatenate((A, b), axis=1) #exindem matricea cu coloana b
    aux = np.zeros((1, lin + 1))
    print(Aext)
    for k in range(lin - 1):
        # valorile de comparare
        max = -1
        p = -1
        # coordonatele elem max din matr pe care le intrschimb in functie de linie
        for i in range(k, lin):
            if abs(Aext[i][k]) > max:
                max = abs(Aext[i][k])
                p = i
        if p == -1:
            print("Sistemul nu admite solutie unica")
            return None, None
        # procesul de interschimbare
        if p != k:
            aux[:] = Aext[k]
            Aext[k] = Aext[p]
            Aext[p] = aux

        # transform elementele sub pivot([k][k]) in 0 folosind transformari elementare pe toate liniile sub linia pivot
        for l in range(k + 1, lin):
            Aext[l] = Aext[l] - (Aext[l][k] / Aext[k][k]) * Aext[k]
        print(Aext)
    if abs(Aext[lin - 1][lin - 1]) <= tol:
        print("Sistemul nu admite solutie unica")
        return None, None
    # stocam x ul calculat
    x = SubsDesc(Aext[:, 0:lin], Aext[:, lin], tol)
    return x, Aext


A = np.array([[1, 5, 3], [5, 13, 12], [1, 2, 1]])
b = np.array([[25], [77], [8]])
tol = 1e-16
xsol, Aext = GaussPivPart(A, b, tol)
print(xsol)
