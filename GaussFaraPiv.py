import numpy as np


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

    for i in range(lin):
        for j in range(i):
            if abs(A[i][j]) > tol:
                print('Matricea nu este superior triunghiulara, dati o alta matrice')
                return None

    for k in range(lin):
        if (abs(A[k][k]) < tol):
            print('Sistemul nu admite solutie unica')
            return None

    x = np.zeros((lin, 1))
    x[lin - 1] = b[lin - 1] / A[lin - 1][lin - 1]
    k = lin - 2
    while (k > - 1):
        sum = 0
        for j in range(k + 1, lin):
            sum += A[k][j] * x[j]
        x[k] = (b[k] - sum) / A[k][k]
        k -= 1
    return x


def GaussFP(A, b, tol):
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

    Aext = np.concatenate((A, b), axis=1)
    aux = np.zeros((1, lin + 1))
    for k in range(lin - 1):
        p = -1
        for i in range(k, lin):
            if abs(Aext[i][k]) > tol:
                p = i
                break
        if p == -1:
            print("Sistemul nu admite solutie unica")
            return None, None

        if p != k:
            aux[:] = Aext[k]
            Aext[k] = Aext[p]
            Aext[p] = aux

        # transform elementele sub pivot([k][k]) in 0 folosind transformari elementare pe toate liniile sub linia pivot
        for l in range(k + 1, lin):
            Aext[l] = Aext[l] - (Aext[l][k] / Aext[k][k]) * Aext[k]

    if abs(Aext[lin - 1][lin - 1]) <= tol:
        print("Sistemul nu admite solutie unica")
        return None, None

    x = SubsDesc(Aext[:, 0:lin], Aext[:, lin], tol)
    return x, Aext

# rezolvare ex 3
# A = np.array([[0, 1, 1], [2, 1, 5], [4, 2, 1]])
# b = np.array([[3], [5], [1]])
# tol = 1e-16
# xsol, Aext = GaussFP(A, b, tol)
# print(xsol)

# print("---Verificare---")
# print(A@xsol)
# print(b)

# #sistem incompatibil, nu admite solutii
# A1 = np.array([[0, 1, -2], [1, -1, 1], [1, 0, -1]])
# b1 = np.array([[4], [6], [2]])
# xsol1, Aext1 = GaussFP(A1, b1, tol)
# print(xsol1)

# print("---Verificare 2---")
# print(A1@xsol1)
# print(b1)