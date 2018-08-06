import sys
sys.path += ['elliptic-curves-finite-fields']

from finitefield.finitefield import FiniteField
from edwards import *

# the order of jubjub base field
q = 0x73eda753299d7d483339d80809a1d80553bda402fffe5bfeffffffff00000001
Fq = FiniteField(q, 1)

dd = -(Fq(10240)/Fq(10241))
curve = TwistedEdwardsCurve(Fq(-1), dd)