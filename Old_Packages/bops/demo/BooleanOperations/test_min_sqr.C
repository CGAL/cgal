//  -*- Mode: c++ -*-
// ============================================================================
// 
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-1.0 $
// release_date  : $CGAL_Date: 1998/09/12 $
//
// file          : demo/BooleanOperations/test_min_sqr.C
// source        : demo/BooleanOperations/test_min_sqr.C
// revision      : $Revision$
// revision_date : $Date$
// author(s)     :                 Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>
//
// coordinator   : RISC Linz
//  (Wolfgang Freiseisen <wfreisei@risc.uni-linz.ac.at>)
//
// 
// ============================================================================

#include <cstdlib>
#include <cmath>
#include <CGAL/min_sqr_distance_traits.h>
#include <CGAL/min_sqr_distance.h>
#include <iostream>
#include <ctime>
#include <sys/timeb.h>
//extern int ftime(struct timeb *__tp);

struct Time {
  typedef unsigned long  seconds;
  typedef unsigned short milliseconds;
  struct timeb _t;

  milliseconds msec() const { return _t.millitm; }
  seconds sec() const { return _t.time; }
  Time() { set(); }
  Time( seconds s, milliseconds m) {
    _t.time= s;
    _t.millitm= m%1000;
  }

  Time( milliseconds m) {
    _t.time= m/1000;
    _t.millitm= m-sec();
  }
  Time(const Time& t) { _t= t._t; }
  Time& operator=(const Time& t) { _t= t._t; return *this; }

  Time operator-( const Time& t) const {
     seconds s= sec() - t.sec() - (msec() < t.msec());
     return Time( s, msec()-t.msec() );
  }
  void set() { ftime(&_t); }
  unsigned long get_msecs() const { return (unsigned long)(sec()*1000 + msec); }
};


ostream& operator<<(ostream& o, const Time& t) {
  o << (float)t.sec() + ((float)t.msec())/1000.;
  return o;
}


typedef double TestNum;  
double Random(double) { return drand48(); }
void Seed(long seed, double) { srand48(seed); }
typedef Cartesian<TestNum>  R_type;
typedef min_sqr_distance_traits<R_type> Traits;
typedef Traits::Point Point;

#include <ctime>

int main(int argc, char *argv[]) {

  list<Point> L;
  int i, n= 24;
  if(argc > 1) n= atoi(argv[1]);
  const double rval= 1000;
  Seed(time(NULL), TestNum(0));
  for(i= 0; i < n; i++)
    L.push_back( Point(Random(TestNum(0))*rval,Random(TestNum(0))*rval) );

  Traits tr;

  cout << "n= " << n << endl;
  //long int t0= time(NULL);
  Time t0;
  t0.set();
  double d1= minimal_square_distance(L.begin(), L.end(), tr);
  //long int t1= time(NULL)-t0;
  Time t1;
  t1= t1 - t0;
  cout << "O(n log n): " << d1 << " (" << t1 << " sec)" << endl;

  t0.set();
  d1= minimal_square_distance(L.begin(), L.end(), tr);
  //long int t1= time(NULL)-t0;
  t1.set();
  t1= t1 - t0;
  cout << "O(n log n): " << d1 << " (" << t1 << " sec)" << endl;


  t0.set();
  d1= minimal_square_distance(L.begin(), L.end(), tr);
  //long int t1= time(NULL)-t0;
  t1.set();
  t1= t1 - t0;
  cout << "O(n log n): " << d1 << " (" << t1 << " sec)" << endl;

  t0.set();
  double d2= minimal_square_distance2(L.begin(), L.end(), tr);
  Time t2;
  //long int t2= time(NULL)-t0;
  t2= t2 - t0;
  cout << "O(n^2): " << d2 << " (" << t2 << " sec)" << endl;

  return 0;
}
