from finitefield.finitefield import FiniteField
from edwards import *
from edwards_proj import *
from edwards_ext import *

import time

# the order of jubjub base field
q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
Fq = FiniteField(q, 1)

dd = -(Fq(10240)/Fq(10241))

curve = TwistedEdwardsCurve(Fq(-1), dd)
curve_proj = ProjectiveEdwards(Fq(-1), dd)
curve_ext = ExtendedEdwards(Fq(-1), dd)

curves = [curve, curve_proj, curve_ext]



def test(curve):
	p = Point(curve, Fq(0x18ea85ca00cb9d895cb7b8669baa263fd270848f90ebefabe95b38300e80bde1), Fq(0x255fa75b6ef4d4e1349876df94ca8c9c3ec97778f89c0c3b2e4ccf25fdf9f7c1))
	q = Point(curve, Fq(0x1624451837683b2c4d2694173df71c9174ffcc613788eef3a9c7a7d0011476fa), Fq(0x6f76dbfd7c62860d59f5937fa66d0571158ff68f28ccd83a4cd41b9918ee8fe2))
	t0 = time.time()
	for i in range(50000):
		p = p + q
	t1 = time.time()
	print curve, "time: ", t1 - t0

for j in range(10):
	print "="*25
	for i in curves:
		test(i)


