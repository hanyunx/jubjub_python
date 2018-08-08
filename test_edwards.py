import sys
sys.path += ['elliptic-curves-finite-fields']

from finitefield.finitefield import FiniteField
from edwards import *

# the order of jubjub base field
q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
Fq = FiniteField(q, 1)

dd = -(Fq(10240)/Fq(10241))
curve = TwistedEdwardsCurve(Fq(-1), dd)

# print dd
# print curve
p1 = Point(curve, 0,1)
# print curve.testPoint(p1.x, p1.y)
p2 = Point(curve, 5, 6846412461894745224441235558443359243034138132682534265960483512729196124138)
# print curve.testPoint(p2.x, p2.y)
# print p1 + p2
# print p2.double().double() + p2
# print 100 * p2
# print p2 - p2
print -p2
print p1 - p2

print 3 * p2 + 6 * p2
print 9 * p2