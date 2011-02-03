#!/usr/bin/python

import subprocess
import re

MAXIMA_FILE = 'equas.mac'

fd = open(MAXIMA_FILE, 'r')

comment = 0
start_comment = 0
end_comment = 0

for line in fd:
  
  start_comment = start_comment + len(re.findall('/\*', line))
  end_comment = end_comment + len(re.findall('\*/', line))

  comment = start_comment - end_comment

  if comment < 0:
      print("Comment syntax error in %s", MAXIMA_FILE)
  elif end_comment:
      line = re.sub('/\*.*\*/', '', line)
  elif comment == 0:
    print(line)

fd.close()
