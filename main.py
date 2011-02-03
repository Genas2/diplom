#!/usr/bin/python

import subprocess

MAXIMA_FILE = 'equas.mac'

fd = open(MAXIMA_FILE, 'r')

comment = False
start_comment = -1
end_comment = -1

for line in fd:
  
  start_comment = line.find('/*')
  end_comment = line.find('*/')

  if start_comment != -1:
    #comment = False
    start_comment = line.find('/*')

    if find('/*') != -1:
      comment = True


  if comment != True:
    print(line)

fd.close()
