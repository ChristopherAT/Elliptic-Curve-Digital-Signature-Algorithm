from scipy import ceil,sqrt,log2
from sympy.ntheory import isprime,primefactors,factorint,nextprime,prevprime
from random import randint
from hashlib import sha1

def PowerMod(x,n,p):
    '''
    Algoritma 2.... (Square-and-Multiply)
    Masukan  : bil bulat x,n,p
    Keluaran : x^n mod p
    '''
    biner=bin(int(n))[2:]
    y=1
    lb=len(biner)-1
    exp=x
    while lb>=0:
        if int(biner[lb])==1:
            y=y*(exp) % p
        exp=(exp**2) %p
        lb-=1
    return y

def Legendre(a,p):
    '''
    Teorema 3.... (Kriteria Euler)
    Masukan  : bil bulat a dan bil prima p
    Keluaran : simbol legendre (a/p)
    '''
    y=PowerMod(a,(p-1)//2,p)
    if y>1:
        return y-p
    else:
        return y

def InverseMod(a,p):
    '''
    Algoritma 2.... (Perluasan Euclidean)
    Masukan  : bil bulat a dan bil prima p
    Keluaran : invers dari a modulo p
    '''
    a=a%p
    u=a
    v=p
    x1=1
    x2=0
    while u>1:
        q=v//u
        r=v-q*u
        x=x2-q*x1

        v=u
        u=r
        x2=x1
        x1=x
    return int(x1)%p

def SqrtMod(a,p):
    '''
    Algoritma 3.... (Akar kuadrat modulo p)
    Masukan  : bil bulat a dan bil prima p
    Keluaran : bil bulat x (x^2=a)
    '''
    leg=Legendre(a,p)
    if leg==-1:
        raise Exception('Tidak ada akar kuadrat modulo ',p,' dari ',a)
    elif leg==0:
        return 0
    elif p%4==3:
        return PowerMod(a,(p+1)//4,p)
    elif p%8==5 and PowerMod(a,(p-1)//4,p)==1:
        return PowerMod(a,(p+3)//8,p)
    elif p%8==5 and PowerMod(a,(p-1)//4,p)==p-1:
        return (2*a*PowerMod(4*a%p,(p-5)//8,p))%p
    else:
        return SqrtModTS(a,p)

def SqrtModTS(a,p):
    '''
    Algoritma 3.... (Tonelli-Shank)
    Masukan  : bil bulat a dan bil prima p
    Keluaran : bil bulat x (x^2=a)
    '''
    q=p-1
    e=0
    while q%2==0:
        e=e+1
        q=q//2

    n=randint(1,p-1)
    while PowerMod(n,((p-1)//2),p)==1:
        n=randint(1,p-1)
    z=PowerMod(n,q,p)

    y=z
    r=e
    x=PowerMod(a,(q-1)//2,p)
    b=(a*x*x)%p
    x=(a*x)%p

    while b!=1:
        m=1
        while PowerMod(b,2**m,p)!=1:
            m=m+1

        if m==r:
            raise Exception('a is not quadratic residue modulo p')
            return -1

        t=PowerMod(y,2**(r-m-1),p)
        y=(t**2)%p
        r=m%p
        x=(x*t)%p
        b=(b*y)%p

    return x

def BabyStepGiantStep(curve,point):
    '''
    Algoritma 3.... (Baby Step Giant Step)
    Masukan  : Kurva eliptik E dan titik P
    Keluaran : Periode titik P
    '''
    P=point
    p=curve.p
    s=int(ceil(p**0.25))
    List_Baby=['Inf']
    Baby=P
    for i in range(1,s+1):
        List_Baby.append(Baby.x)
        Baby=Baby+P

    Q=(2*s+1)*P
    R=(p+1)*P
    j=0
    jQ=Infty(curve)
    while True:
        Giant=R+jQ
        if Giant.x in List_Baby:
            i=List_Baby.index(Giant.x)
            Baby=i*P
            if Giant.x==Baby.x:
                return (p+1)+(2*s+1)*j-i
            else:
                return (p+1)+(2*s+1)*j+i
        Giant=R-jQ
        if Giant.x in List_Baby:
            i=List_Baby.index(Giant.x)
            Baby=i*P
            if Giant.x==Baby.x:
                return (p+1)-(2*s+1)*j-i
            else:
                return (p+1)-(2*s+1)*j+i
        jQ=jQ+Q
        j+=1

def OrderTitik(Periode,P):
	Inf=P-P
	faktor=primefactors(Periode)
	for i in faktor:
		if (Periode//i)*P==Inf:
			return OrderTitik(Periode//i,P)
	return Periode

def OrderGrupNaive(curve):
    '''
    Algoritma 3.... (Algoritma Naive)
    Masukan  : Kurva eliptik E
    Keluaran : Order grup E(Z_p)
    '''
    order=curve.p+1
    for x in range (curve.p):
        y=(x**3+curve.a*x+curve.b)%curve.p
        order=order+Legendre(y,curve.p)
    return order

def OrderGrupShank(kurva,P=None,coba=None):
    '''
    Algoritma 3.... (Algoritma Shank)
    Masukan  : Kurva eliptik E
    Keluaran : Order grup E(Z_p)
    '''
    p=kurva.p
    if coba==None:
        coba=int(log2(float(p)))+5
    if P==None:
        P=GeneratorTitik(kurva)

    M=BabyStepGiantStep(kurva,P)
    M=OrderTitik(M,P)

    List_M=[]
    N=int(ceil(((p+1-2*sqrt(p))/M)))*M
    while N<=p+1+2*sqrt(p):
        List_M.append(N)
        N=N+M

    counter=0
    while len(List_M)!=1 and counter<coba:
        Q=GeneratorTitik(kurva)
        for Mi in List_M:
            MiQ=Mi*Q
            if MiQ!=Infty(kurva):
                List_M.remove(Mi)
        counter+=1

    if counter==coba:
        print('Error: order not found')
        return 0
    return List_M[0]


def hash(s):
    ''' Output nilai hash dari s. Jika input bukan berupa bytes,
        akan di-encode terlebih dahulu dengan ascii.
    '''
    if type(s) is bytes:
        return sha1(s).hexdigest()
    else:
        return sha1(s.encode('ascii')).hexdigest()

def PollardRho(n,P,Q):
    '''
    Algoritma Pollard-Rho
    Masukan : titik P,Q dan n = order P
    Keluaran: m=log_P(Q)
    '''
    if not isprime(n):
        print('Order titik P tidak prima!')
    def f(Z,a,b):
        if Z==Inf:
            return (P+Z,(a+1)%n,b)
        elif Z.x%3==1:
            return (P+Z,(a+1)%n,b)
        elif Z.x%3==0:
            return (Z+Z,(2*a)%n,(2*b)%n)
        elif Z.x%3==2:
            return (Q+Z,a,(b+1)%n)
    Inf=P-P
    (Zi,ai,bi)=f(Q,0,1)
    (Z2i,a2i,b2i)=f(Zi,ai,bi)
    while True:
        (Zi,ai,bi)=f(Zi,ai,bi)
        (Z2i,a2i,b2i)=f(Z2i,a2i,b2i)
        (Z2i,a2i,b2i)=f(Z2i,a2i,b2i)
        if Zi==Z2i:
            if b2i-bi==0:
                return 0
            else:
                m=((ai-a2i)*InverseMod(b2i-bi,n)) % n
                return m

def GeneratorKunci(DP):
    privat=randint(1,DP.order)
    publik=privat*DP.generator
    return (privat,publik)

def TandaTangan(x,privat,DP):
    e=int(hash(x),16)
    k=randint(1,DP.order-1)
    kA=k*DP.generator
    u,v=kA.x,kA.y
    r=u%DP.order
    s=InverseMod(k,DP.order)*(e+privat*r) % DP.order
    if r==0 or s==0:
        return TandaTangan(x,DP.generator*privat,DP)
    else:
        return (r,s)

def Verifikasi(x,ttd,publik,DP):
    r=ttd[0]
    s=ttd[1]
    sha1x=hash(x)
    w=InverseMod(s,DP.order)
    i=w*int(sha1x,16) % DP.order
    j=w*r % DP.order
    uv=i*DP.generator+j*publik
    if uv.x % DP.order == r:
        return True
    else:
        return False

class KurvaEliptik(object):
   def __init__(self, a, b, p):
      # kurva: y^2=x^3+ax+b
      self.a = a
      self.b = b
      self.p = p

      self.diskriminan = -16*(4 * a*a*a + 27 * b * b) % p
      if not self.isSmooth():
         raise Exception("Kurva eliptik %s tidak smooth!" % self)

   def isSmooth(self):
      return self.diskriminan != 0

   def testPoint(self, x, y):
      return (y*y) % self.p == (x*x*x + self.a * x + self.b) % self.p

   def __str__(self):
      return 'y^2 = x^3 + %sx + %s atas Z_%s' % (self.a, self.b, self.p)

   def __repr__(self):
      return str(self)

   def __eq__(self, other):
      return (self.a, self.b, self.p) == (other.a, other.b, self.p)


class Titik(object):
   def __init__(self, kurva, x, y):
      self.kurva = kurva
      self.x = x % kurva.p
      self.y = y % kurva.p

      if not kurva.testPoint(self.x,self.y):
         raise Exception("Titik %s tidak berada pada kurva eliptik %s!" % (self, kurva))

   def __str__(self):
      return "(%r, %r)" % (self.x, self.y)

   def __repr__(self):
      return str(self)

   def __neg__(self):
      return Titik(self.kurva, self.x, self.kurva.p-self.y)

   def __add__(self, Q):
      if self.kurva != Q.kurva:
         raise Exception("Tidak bisa menjumlahkan dua titik dengan kurva eliptik yang berbeda!")
      if isinstance(Q, Infty):
         return self

      x_1, y_1, x_2, y_2 = self.x, self.y, Q.x, Q.y
      a,p=self.kurva.a, self.kurva.p

      if (x_1, y_1) == (x_2, y_2):
         if y_1 == 0:
            return Infty(self.kurva)

         m = ((3 * x_1 * x_1 + a) * InverseMod(2 * y_1,p)) % p
      else:
         if x_1 == x_2:
            return Infty(self.kurva)

         m = ((y_2 - y_1) * InverseMod(x_2 - x_1,p)) % p

      x_3 = (m*m - x_1 - x_2) % p
      y_3 = (m*(x_1 - x_3) - y_1) % p

      return Titik(self.kurva, x_3, y_3)


   def __sub__(self, Q):
      return self + -Q

   def __mul__(self, n):
      if not (isinstance(n, int) or isinstance(n, long)):
         raise Exception("Perkalian harus dengan bilangan bulat!")
      else:
         if n < 0:
             return -self * -n
         if n == 0:
             return Infty(self.kurva)
         else:
             Q = self
             R = self if n & 1 == 1 else Infty(self.kurva)
             i = 2
             while i <= n:
                 Q = Q + Q
                 if n & i == i:
                     R = Q + R
                 i = i << 1
             return R


   def __rmul__(self, n):
      return self * n

   def __list__(self):
      return [self.x, self.y]

   def __eq__(self, other):
      if type(other) is Infty:
         return False

      return (self.x, self.y) == (other.x, other.y)

   def __ne__(self, other):
      return not self == other

   def __getitem__(self, index):
      return [self.x, self.y][index]

def GeneratorTitik(kurva):
        x=randint(0,kurva.p-1)
        while True:
            x2=(x*x*x+kurva.a*x+kurva.b)%kurva.p
            if Legendre(x2,kurva.p)==1:
                y=SqrtMod(x2,kurva.p)
                break
            else:
                x=(x+1)%kurva.p
        return Titik(kurva,x,y)


class Infty(Titik):
   def __init__(self, kurva):
      self.kurva = kurva
      self.x='Inf'
      self.y='Inf'

   def __neg__(self):
      return self

   def __str__(self):
      return "Infty"

   def __add__(self, Q):
      if self.kurva != Q.kurva:
         raise Exception("Tidak bisa menjumlahkan dua titik dengan kurva eliptik yang berbeda!")
      return Q

   def __mul__(self, n):
      if not (isinstance(n, int) or isinstance(n, long)):
         raise Exception("Perkalian harus dengan bilangan bulat!")
      else:
         return self

   def __eq__(self, other):
      return type(other) is Infty

class DomainParameter(object):
    def __init__(self, kurva, P, order):
      self.kurva=kurva
      self.generator=P
      self.order=order

    def __str__(self):
      return 'Kurva eliptik   : %s \nTitik generator : %s\nOrder grup      : %s' % (self.kurva, self.generator, self.order)

    def __repr__(self):
      return str(self)

    def __eq__(self, other):
      return (self.kurva, self.generator) == (other.kurva, other.generator)

class DomainParameterNIST(DomainParameter):
    def __init__(self):
        p=2**192-2**64-1
        a=-3
        b=0x64210519E59C80E70FA7E9AB72243049FEB8DEECC146B9B1
        n=0xFFFFFFFFFFFFFFFFFFFFFFFF99DEF836146BC9B1B4D22831
        x=0x188DA80EB03090F67CBF20EB43A18800F4FF0AFD82FF1012
        y=0x07192B95FFC8DA78631011ED6B24CDD573F977A11E794811

        self.kurva=KurvaEliptik(a,b,p)
        self.generator=Titik(self.kurva,x,y)
        self.order=n