#!/usr/bin/python

def legendre_calculate(n, x):
  """ Calculates n-th Legendre polynomial """

  result = 1
  p = [1, x]

  s = ''

  if n == 0:
   result = p[0]
  elif n == 1:
   result = p[1]
  else:
   for i in range(1,n+1):
    p.append(((2 * i - 1) * x * p[-1] - (i-1) * p[-2]) / i)

  return p[-1]

#print(legendre_calculate(2,1))
