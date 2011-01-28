#!/usr/bin/python

import math

def diff(func, x, l, exponent=0, prec=0.01):
  if exponent == 0:
   return func(x, l)
  elif exponent < 0:
   return -1
  elif exponent == 1:
   return (func(x + prec, l) - func(x - prec, l))/(2 * prec)

  while exponent > 1:
   expontent = exponent - 1
   result = (diff(func, x + prec, l, exponent, prec) - diff(func, x - prec, l, exponent, prec))/(2 * prec) 

func = lambda x, l: math.pow((math.cos(x) * math.cos(x) - 1), l)

print(diff(func, math.pi / 4, 2, 1))

