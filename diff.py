#!/usr/bin/python

import math

def diff(func, x, l, prec=0.01):
  return (func(x + prec, l) - func(x - prec, l))/(2 * prec)

func = lambda x, l: math.pow((math.cos(x) * math.cos(x) - 1), l)

print(diff(func, math.pi/2, 1))

