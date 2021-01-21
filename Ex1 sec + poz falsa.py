import numpy as np
import matplotlib.pyplot as plt

e=pow(10,-5) #epsilon
f = lambda x: x**3 - x**2 - 10*x + 10 #functia mea

def MetSecantei(f, x0, x1, e=pow(10,-5)):
    nr_iteratii = 0 #numarul de pasi

    while abs(x1 - x0) / abs(x0) >= e:
        nr_iteratii = nr_iteratii + 1
        #trecem la urm valoare in functie de x1 si x0
        x_urm = x1 - (f(x1)*(x1-x0))/(f(x1)-f(x0))
        x0 = x1
        x1 = x_urm

        if (x_urm < a or x_urm > b):
            break

    return x_urm, nr_iteratii
#capetele intervalelor
a = float(input("a = "))
b = float(input("b = "))

#cele 3 radacini
rad1, nr_it1 = MetSecantei(f,  0.1, 2)
rad2, nr_it2 = MetSecantei(f,  -1, 5)
rad3, nr_it3 = MetSecantei(f,  2, 4)
print("Solutia 1:  %0.5f , Numarul de pasi:  %d " %(rad1, nr_it1))
print("Solutia 2:  %0.5f , Numarul de pasi:  %d " %(rad2, nr_it2))
print("Solutia 3:  %0.5f , Numarul de pasi:  %d " %(rad3, nr_it3))

#metoda pozitie false
def MetPozFalse(a,b,e, f):
    nr_iteratii= 0 #numarul de pasi
    cond = True
    while cond:
        nr_iteratii = nr_iteratii + 1
        #calculam radacina in functie de capetele intervalelor
        c = a - (b-a) * f(a)/( f(b) - f(a) )
        #verificam care capat indeplineste conditia f (a)f (b) < 0 pt a determina
        # noul capat "c" al intervalului || verificam daca sunt de semne dif
        if f(a) * f(c) < 0:
            b = c
        else:
            a = c
        #conditia este ca functia de c sa fie mai mare ca epsilon
        cond = abs(f(c)) >= e
    return c,nr_iteratii




if f(a) * f(b) > 0.0:
   print('Functia nu ideplineste conditia')
else:
    c,nr_iteratii = MetPozFalse(a,b,e,f)
    print('Radacina necesara este: %0.5f cu numarul de iteratii: %d' % (c, nr_iteratii))
plt.scatter([rad1,rad2,rad3],[0,0,0],c='m')
x = np.linspace(-10, 10)
y =  x**3 - x**2 - 10*x + 10
plt.plot(x,y,'m')
plt.show()