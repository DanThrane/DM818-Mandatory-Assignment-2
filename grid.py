#!/usr/bin/python

for x in range(0, 71):
  for y in range(0, 71):
    print format(y + x * 71, '04'),
  print ""
