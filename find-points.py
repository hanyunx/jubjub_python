import sys
sys.path += ['elliptic-curves-finite-fields']

from edwards import *
from finitefield.finitefield import FiniteField
from fractions import Fraction as frac

import itertools

def findPoints(curve, field):
   print('Finding all points over %s' % (curve))
   print('The ideal generator is %s' % (field.idealGenerator))

   degree = field.idealGenerator.degree()
   subfield = field.primeSubfield
   xs = [field(x) for x in itertools.product(range(subfield.p), repeat=degree)]
   ys = [field(x) for x in itertools.product(range(subfield.p), repeat=degree)]

   points = [Point(curve, x, y) for x in xs for y in ys if curve.testPoint(x,y)]
   return points


q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
Fq = FiniteField(q, 1)

dd = -(Fq(10240)/Fq(10241))
curve = TwistedEdwardsCurve(Fq(-1), dd)
points = findPoints(curve, Fq)

for point in points:
   print(point)
