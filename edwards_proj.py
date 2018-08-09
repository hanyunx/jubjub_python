from finitefield.finitefield import FiniteField

q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
Fq = FiniteField(q, 1)

# Twisted Edwards Curve
class ProjectiveEdwards(object):
   def __init__(self, a, d):
      self.a = a
      self.d = d

      self.disc = a * d * (a - d) * (a - d) * (a - d) * (a - d)
      self.j = 16 * (a * a + 14 * a * d + d * d) * (a * a + 14 * a * d + d * d) * \
               (a * a + 14 * a * d + d * d) / self.disc
      if not self.isSmooth():
         raise Exception("The curve %s is not smooth!" % self)


   def isSmooth(self):
      return self.disc != 0


   def testPoint(self, x, y):
      return self.a * x * x + y*y == 1 + self.d * x * x * y * y


   def __str__(self):
      return '%sx^2 + y^2 = 1 + %sx^2y^2' % (self.a, self.d)


   def __repr__(self):
      return str(self)


   def __eq__(self, other):
      return (self.a, self.d) == (other.a, other.d)



class Point(object):
   def __init__(self, curve, x, y):
      self.curve = curve # the curve containing this point
      self.x = x
      self.y = y

      if not curve.testPoint(x,y):
         raise Exception("The point %s is not on the given curve %s!" % (self, curve))


   def __str__(self):
      return "(%r, %r)" % (self.x, self.y)


   def __repr__(self):
      return str(self)


   def __neg__(self):
      return Point(self.curve, -self.x, self.y)


   def __add__(self, Q):
      # Assumptions: Z1=1 and Z2=1.
      # Cost: 6M + 1S + 1*a + 1*d + 8add.
      # source 2008 Bernstein--Birkner--Joye--Lange--Peters http://eprint.iacr.org/2008/013 Section 6, plus Z2=1, plus Z1=1, plus standard simplification
      # assume Z1 = 1, Z2 = 1
      # compute C = X1 X2
      #         D = Y1 Y2
      #         E = d C D
      #         X3 = (1-E) ((X1+Y1)(X2+Y2)-C-D)
      #         Y3 = (1+E) (D-a C)
      #         Z3 = 1-E^2

      if self.curve != Q.curve:
         raise Exception("Can't add points on different curves!")
      if isinstance(Q, Ideal):
         return self

      x1, y1, z1 = self.x, self.y, 1
      x2, y2, z2 = Q.x, Q.y, 1

      c = x1 * x2
      d = y1 * y2
      e = self.curve.d * c * d

      x3 = (1 - e) * ((x1 + y1) * (x2 + y2) - c - d)
      y3 = (1 + e) * (d + c)
      z3 = 1 - e * e
      
      return Point(self.curve, x3/z3, y3/z3)


   def double(self):
      # Assumptions: Z1=1.
      # Cost: 2M + 4S + 1*a + 7add + 1*2.
      # https://hyperelliptic.org/EFD/g1p/auto-twisted-projective.html#doubling-mdbl-2008-bbjlp
      # assume Z1 = 1
      # compute B = (X1+Y1)^2
      #         C = X1^2
      #         D = Y1^2
      #         E = a C = -C
      #         F = E + D
      #         X3 = (B-C-D)(F-2)
      #         Y3 = F(E-D)
      #         Z3 = F^2-2F

      x1, y1, z1 = self.x, self.y, 1

      b = (x1 + y1) * (x1 + y1)
      c = x1 * x1
      d = y1 * y1
      e = -c
      f = e + d

      x3 = (b - c - d) * (f - 2)
      y3 = f * (e - d)
      z3 = f * f - 2 * f

      X3 = x3*(Fq(z3).inverse())
      Y3 = y3*(Fq(z3).inverse())

      return Point(self.curve, X3, Y3)
      # return Point(self.curve, x3/z3, y3/z3)
      # return self + self


   def __sub__(self, Q):
      return self + -Q

   def __mul__(self, n):
      if not isinstance(n, int):
         raise Exception("Can't scale a point by something which isn't an int!")

      if n < 0:
         return -self * -n

      if n == 0:
         return Ideal(self.curve)

      Q = self
      R = self if n & 1 == 1 else Ideal(self.curve)

      i = 2
      while i <= n:
         Q += Q

         if n & i == i:
             R += Q

         i = i << 1

      return R


   def __rmul__(self, n):
      return self * n

   def __list__(self):
      return [self.x, self.y]

   def __eq__(self, other):
      if type(other) is Ideal:
         return False

      return self.x, self.y == other.x, other.y

   def __ne__(self, other):
      return not self == other

   def __getitem__(self, index):
      return [self.x, self.y][index]

# TODO?
class Ideal(Point):
   def __init__(self, curve):
      self.curve = curve

   def __neg__(self):
      return self

   def __str__(self):
      return "Ideal"

   def __add__(self, Q):
      if self.curve != Q.curve:
         raise Exception("Can't add points on different curves!")
      return Q

   def __mul__(self, n):
      if not isinstance(n, int):
         raise Exception("Can't scale a point by something which isn't an int!")
      else:
         return self

   def __eq__(self, other):
      return type(other) is Ideal

