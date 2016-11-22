#!/usr/bin/env python3.4

import sys
import argparse

class Point_2:
  def __init__(self, x, y):
    self.x = int(x)
    self.y = int(y)
  def __str__(self):
    return 'Point_2(' + str(self.x) + ', ' + str(self.y) + ')'

class Segment_2:
  def __init__(self, source, target):
    self.source = source
    self.target = target
  def __str__(self):
    return 'Segment_2(' + str(self.source) + ', ' + str(self.target) + ')'

class Sign:
  def __init__(self, sgn):
    assert(sgn in [-1, 0, +1])
    self.s = sgn
  def __str__(self):
    if self.s == -1:
      return 'CGAL::NEGATIVE'
    elif self.s == 0:
      return 'CGAL::ZERO'
    elif self.s == +1:
      return 'CGAL::POSITIVE'
    else:
      assert(False)

def clean(s):
  if s[-1] == ')':
    return clean(s[:-1])
  if s[0] == '[':
    res = s.split(';')
    assert(len(res) == 2)
    assert(res[0][1:] == res[1][:-1])
    return res[0][1:]
  else:
    return s

def print_test(p, q, r, t, sgn):
  print('  test_incircle<Gt>(')
  print('      ', p, end='')
  print(',')
  print('      ', q, end='')
  print(',')
  print('      ', r, end='')
  print(',')
  print('      ', t, end='')
  print(',')
  print('      ', sgn, end='')
  print(');')

def test_sp(s):
  return ((s[-1] == 'p') or (s[-1] == 's')) and len(s) < 4

def count_points(sites):
  count = 0
  for i in [0, 1, 2]:
    if isinstance(s[i], Point_2):
      count = count + 1
  return count

# main

parser = argparse.ArgumentParser(\
    description='Convert vertex conflict return values to test format.')
parser.add_argument('--points')
args = parser.parse_args()
#print(args)
if args.points:
  num_points = int(args.points)
  assert(num_points in [0, 1, 2, 3])

for line in sys.stdin:
  a = line.strip().split()
  #print(a)
  lenarray = len(a)
  if lenarray < 4*3+1:
    continue
  i = 0
  sites_read = 0
  s = []
  while (sites_read < 4):
    #print(sites_read)
    while (i < lenarray) and (not test_sp(a[i])):
      i = i + 1
    assert(a[i][-1] == 'p' or a[i][-1] == 's')
    is_point = a[i][-1] == 'p'
    #print(i)
    if (is_point):
      p = clean(a[i+1])
      q = clean(a[i+2])
      s.append( Point_2(p, q) )
      i = i + 3
    else:
      p1 = clean(a[i+1])
      q1 = clean(a[i+2])
      p2 = clean(a[i+3])
      q2 = clean(a[i+4])
      s.append( \
        Segment_2(Point_2(p1, q1), Point_2(p2, q2)) )
      i = i + 5
    sites_read = sites_read + 1

  assert(a[-1] in ['-1', '0', '+1', '1'])

  if args.points:
    count = count_points(s)
    #print(count)

  if (not args.points) or (count == num_points):
    print_test(s[0], s[1], s[2], s[3], Sign(int(a[-1])))
