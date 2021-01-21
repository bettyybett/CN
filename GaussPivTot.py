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
    # if (lin != col):
    #     col = np.shape(A)[1]
    #     print('Matricea nu este patratica, dati alta matrice')
    #     return None
    # print(A)
    # for i in range(lin):
    #     for j in range(i):
    #         print(A[i][j])
    #         if abs(A[i][j])  > tol:
    #             print('Matricea nu este superior triunghiulara, dati o alta matrice')
    #             return None

    for k in range(lin):
        if (abs(A[k][k]) < tol):
            print('Sistemul nu admite solutie unica')
            return None
    #calculam x ul nostru din Ax=b
    x = np.zeros((lin, 1))
    x[lin - 1] = b[lin - 1] / A[lin - 1][lin - 1] #calculam ultima valoare a lui x
    k = lin - 2 #valoare de la care se calculeaza x ul
    while (k > - 1):
        sum = 0
        for j in range(k + 1, lin):
            sum += A[k][j] * x[j]
        x[k] = (b[k] - sum) / A[k][k]
        k -= 1
    return x
#functia de determ a matricei
def matrice(d,f,c,n):
    A = np.zeros((n,n))
    A[0][0] = d
    A[0][1] = f
    for i in range(1,n-1):
        A[i][i-1] = c
        A[i][i] = d
        A[i][i+1] = f
    A[n-1][n-2] = c
    A[n-1][n-1] = d
    return A



def GaussPivTot(A, b, tol):
    '''
    GaussPivTot rezolva sisteme de ecuatii liniare folosind metoda Gauss cu pivotare totala
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
    index = [x for x in range(lin)]
    x = np.zeros((lin, 1))
    aux = np.zeros((1, lin + 1))
    aux1 = np.zeros((lin + 1, 1))

    for k in range(lin - 1):
        #valorile de comparare
        max = -1
        p = -1
        m = -1
        #coordonatele elem max din matr pe care le intrschimb in functie de coloana/linie
        for i in range(k, lin):
            for j in range(k, col):
                if abs(Aext[i][j]) > max:
                    max = abs(Aext[i][j])
                    p = i
                    m = j

        if p == -1 and m == -1:
            print("Sistemul nu admite solutie unica")
            return None, None
        #procesul de interschimbare
        if p != k:
            aux[:] = Aext[k]
            Aext[k] = Aext[p]
            Aext[p] = aux
        if m != k:
            aux1[:, 0] = Aext[:][k]
            Aext[:][k] = Aext[:][m]
            Aext[:][m] = aux1[:, 0]
            aux2 = index[k]
            index[k] = index[m]
            index[m] = aux2

        # transform elementele sub pivot([k][k]) in 0 folosind transformari elementare pe toate liniile sub linia pivot
        for l in range(k + 1, lin):
            Aext[l] = Aext[l] - (Aext[l][k] / Aext[k][k]) * Aext[k]

    if abs(Aext[lin - 1][lin - 1]) <= tol:
        print("Sistemul nu admite solutie unica")
        return None, None
    #stocam x ul calculat
    x = SubsDesc(Aext[:, 0:lin], Aext[:, lin], tol)
    #print(y)
    #stocam solutiile din y in x
    # for i in range(len(y)):
    #     x[index[i]] = y[i]
    return x, Aext


A = np.array(matrice(21,-7,-2,20))#np.array([[0, 1, 1], [2, 1, 5], [4, 2, 1]])
b = np.array([[2],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[1],[2]])
tol = 1e-16

xsol, Aext = GaussPivTot(A, b, tol)
print(xsol)