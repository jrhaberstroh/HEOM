# Parameters using spectroscopic units:
#  2 pi c = 1
#  hc     = 1
import math as m

hbar= 1./((2. * m.pi) * (2. * m.pi))  #It's just true.
hbarsq= hbar * hbar
invhbar = 1/hbar
invhbarsq = 1/hbarsq


damp=1./50. * 1./5309. #fs -> cm
reorg=35.              #cm-1
kt=300. * .695         #k -> cm-1
