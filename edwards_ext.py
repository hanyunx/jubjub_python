from finitefield.finitefield import FiniteField

q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
Fq = FiniteField(q, 1)

# Twisted Edwards Curve
class ExtendedEdwards(object):
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
      # https://hyperelliptic.org/EFD/g1p/auto-twisted-extended.html
      # Assumptions: Z1=1 and Z2=1.
      # Cost: 7M + 1S + 1*a + 1*d + 8add.
      # Cost: 7M + 1S + 1*a + 7add dependent upon the first point.
      # Strongly unified.

      if self.curve != Q.curve:
         raise Exception("Can't add points on different curves!")
      if isinstance(Q, Ideal):
         return self

      x1, y1, t1, z1 = self.x, self.y, self.x * self.y, 1
      x2, y2, t2, z2 = Q.x, Q.y, Q.x * Q.y, 1

      a = x1 * x2
      b = y1 * y2
      c = self.curve.d * t1 * t2
      d = z1 * z2
      h = b + a
      e = (x1 + y1) * (x2 + y2) - h
      f = 1 - c
      g = 1 + c

      x3 = e * f
      y3 = g * h
      t3 = e * h
      z3 = f * g
      
      return Point(self.curve, x3/z3, y3/z3)

   def double(self):
      # See "Twisted Edwards Curves Revisited" Section 3.3
      # http://hyperelliptic.org/EFD/g1p/auto-twisted-extended.html#doubling-dbl-2008-hwcd
      # Cost: 4M + 4S + 1*a + 6add + 1*2.
      # assume: Z = 1
      
      a = self.x * self.x
      b = self.y * self.y
      c = 2
      d = -a
      e = (self.x + self.y) * (self.x + self.y) - a - b
      g = d + b
      f = g - c
      h = d - b

      x3 = e * f
      y3 = g * h
      t3 = e * h
      z3 = f * g

      X3 = x3*(Fq(z3).inverse())
      Y3 = y3*(Fq(z3).inverse())

      return Point(self.curve, X3, Y3)
      # return Point(self.curve, x3/z3, y3/z3)


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

