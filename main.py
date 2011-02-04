#!/usr/bin/python

import subprocess
import re

MAXIMA_FILE = 'equas.mac'

def read_maxima_file(maxima_file):

  equations = {}

  fd = open(maxima_file, 'r')
  
  comment = 0
  start_comment = 0
  end_comment = 0
  
  for line in fd:
    
    start_comment = len(re.findall('/\*', line))
    end_comment = len(re.findall('\*/', line))
  
    comment = comment + start_comment - end_comment
  
    if comment < 0:
      print("Comment syntax error in " + maxima_file)
    elif start_comment and end_comment:
      line = re.sub('/\*[^\*/]*\*/', '', line)
    elif start_comment:
      line = re.sub('/\*[^/\*]*', '', line)
    elif end_comment:
      line = re.sub('[^\*/]*\*/', '', line)
    elif comment != 0:
      continue
  
    if line and line != '\n':
      eq_name = line[0:line.find('(')] 
      equations[eq_name] = line

  fd.close()
  return equations

print(read_maxima_file(MAXIMA_FILE))
