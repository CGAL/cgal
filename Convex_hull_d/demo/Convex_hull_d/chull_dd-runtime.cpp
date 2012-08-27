// Copyright (c) 2001, 2002, 2004  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#include <CGAL/basic.h>
#include <CGAL/Homogeneous_d.h>
#include <CGAL/Cartesian_d.h>
#include <CGAL/Convex_hull_d.h>
#include <CGAL/random_selection.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>


#ifdef CGAL_USE_GMP
#include <CGAL/Gmpz.h>
typedef CGAL::Homogeneous_d<CGAL::Gmpz>   GKernel;
typedef CGAL::Convex_hull_d<GKernel> GConvex_hull_d;
#endif

#ifdef CGAL_USE_LEDA
#include <CGAL/leda_integer.h>
#include <CGAL/leda_real.h>
typedef CGAL::Homogeneous_d<leda_integer> LHKernel;
typedef CGAL::Convex_hull_d<LHKernel> LHConvex_hull_d;
typedef CGAL::Cartesian_d<leda_real> LCKernel;
typedef CGAL::Convex_hull_d<LCKernel> LCConvex_hull_d;
#endif

#include <CGAL/double.h>
typedef CGAL::Homogeneous_d<double>   DHKernel;
typedef CGAL::Convex_hull_d<DHKernel> DHConvex_hull_d;
typedef CGAL::Cartesian_d<double>     DCKernel;
typedef CGAL::Convex_hull_d<DCKernel> DCConvex_hull_d;
typedef int* p_int;
typedef p_int* pp_int;

static std::ofstream* p_table_file;
CGAL::Timer timer;

void random_d_tuple_in_range(p_int V, int d, int l, int h)
{
  for(int i = 0; i<d; ++i)
    V[i] = CGAL::default_random.get_int(l,h);
}

void random_d_tuples_in_range(pp_int V, int n, int d, int l, int h)
{ for(int i = 0; i<n; ++i)
    random_d_tuple_in_range(V[i],d,l,h);
}

void create(pp_int& V, int n, int d)
{
  V = new p_int[n];
  for (int i=0; i<n; ++i)
    V[i] = new int[d];
}

void destroy(pp_int& V, int n)
{ for (int i=0; i<n; ++i) delete [] V[i];
  delete [] V;
  V = 0;
}

void print_to_file(pp_int V, int n, int d, std::string s)
{
  std::ofstream Out(s.c_str());
  Out << d << std::endl;
  Out << n << std::endl;
  for (int i=0; i<n; ++i) {
    for (int j=0; j<d; ++j) Out << V[i][j] << " ";
    Out << std::endl;
  }
  Out.close();
}

void read_from_file(pp_int& V, int& n, int& d, std::string s)
{
  std::ifstream In(s.c_str());
  CGAL_assertion(In.good());
  In >> d >> n;
  create(V,n,d);int i=0,c;
  while ( In >> c ) { V[i/d][i%d]=c; ++i;}
  In.close();
}

template <class R>
void time_insertion_and_check(pp_int V, int n, int d,
  CGAL::Convex_hull_d<R>& C, std::string s, bool check=true)
{
  typedef typename CGAL::Convex_hull_d<R>::chull_has_local_non_convexity
    chull_has_local_non_convexity;
  typedef typename CGAL::Convex_hull_d<R>::chull_has_double_coverage
    chull_has_double_coverage;
  typedef typename CGAL::Convex_hull_d<R>::
    chull_has_center_on_wrong_side_of_hull_facet
    chull_has_center_on_wrong_side_of_hull_facet;

  std::cout << " timing of " << s << std::endl;
  std::vector< CGAL::Point_d<R> > P(n); int i;
  for(i=0; i<n; ++i)
    P[i] = CGAL::Point_d<R>(d,V[i],V[i]+d,1);

  timer.reset(); timer.start(); // float ti = used_time();
  for(i=0; i<n; ++i) {
    C.insert(P[i]);
    if (i%10==0) std::cout << i << " points inserted" << std::endl;
  }
  timer.stop();
  double t = timer.time(); timer.reset(); // float t = used_time(ti);
  (*p_table_file) << s << "\t" << d << " " << n << " "
                  << C.number_of_vertices() << " " << C.number_of_facets()
                  << "\t" << t;
  C.print_statistics();
  std::cout << "used time for inserts  " << t << std::endl;

  C.clear(d);
  timer.start(); // ti = used_time();
  C.initialize(P.begin(),P.end());
  timer.stop(); t = timer.time(); timer.reset();
  // t = used_time(ti);
  C.print_statistics();
  std::cout << "used time for inserts  " << t << std::endl;

  if (check) {
    timer.start();
    std::cout << "entering check" << std::endl;
    try { C.is_valid(true); }
    catch ( chull_has_local_non_convexity )
    { std::cerr << "local non-convexity determined\n"; }
    catch ( chull_has_double_coverage )
    { std::cerr << "double coverage determined\n"; }
    catch ( chull_has_center_on_wrong_side_of_hull_facet )
    { std::cerr << "facet center problem determined\n"; }

    // t = used_time(ti);
    timer.stop(); t = timer.time();
    (*p_table_file) << "\t" << t <<std::endl;
    std::cout<<"used time for sanity check  "<< t <<std::endl<<std::endl;
  } else {
    (*p_table_file) << "\t" << "no"<<std::endl;
    std::cout<<"no check"<<std::endl;
  }
  p_table_file->flush();
}


int main(int argc, char* argv[])
{
  // first param is dimension
  // second param is number of points
  int d = 4;
  int n = 100;
  int range = 100;
  unsigned which(255);

  if (argc > 1 && std::string(argv[1])=="-h") {
    std::cout << "usage: " << argv[0] << " [which] [dim] [#points] [range]\n";
    std::cout << "       which=0 create testset in "<<argv[0]<<".ch\n" ;
    std::cout << "       which=1 double cartesian\n" ;
    std::cout << "       which=2 double homogeneous\n" ;
    std::cout << "       which=4 LEDA integer homogeneous\n" ;
    std::cout << "       which=8 GNU mpz homogeneous\n" ;
    std::cout << "       which=16 LEDA real cartesian\n" ;
    std::exit(1);
  }
  if (argc > 1) which = std::atoi(argv[1]);
  if (argc > 2) d = std::atoi(argv[2]);
  if (argc > 3) n = std::atoi(argv[3]);
  if (argc > 4) range = std::atoi(argv[4]);
  p_table_file = new std::ofstream(
    (std::string(argv[0])+".rts").c_str(), std::ios::app);

  int** V;
  if ( which == 0 ) {
    create(V,n,d);
    random_d_tuples_in_range(V,n,d,-range,range);
    print_to_file(V,n,d,std::string(argv[0])+".ch");
    std::exit(0);
  } else {
    read_from_file(V,n,d,std::string(argv[0])+".ch");
  }

if ( which & 1 ) {
  DCConvex_hull_d DCC(d);
  time_insertion_and_check(V,n,d,DCC,"double cartesian        ");
}
if ( which & 2 ) {
  DHConvex_hull_d DHC(d);
  time_insertion_and_check(V,n,d,DHC,"double homogeneous      ");
}

#ifdef CGAL_USE_LEDA
if ( which & 4 ) {
  LHConvex_hull_d LHC(d);
  time_insertion_and_check(V,n,d,LHC,"LEDA integer homogeneous");
}
#endif
#ifdef CGAL_USE_GMP
if ( which & 8 ) {
  GConvex_hull_d GC(d);
  time_insertion_and_check(V,n,d,GC,"GNU mpz homogeneous     ");
}
#endif
#ifdef CGAL_USE_LEDA
if ( which & 16 )  {
  LCConvex_hull_d LCC(d);
  time_insertion_and_check(V,n,d,LCC,"LEDA real cartesian     ",false);
}
#endif

#ifdef CGAL_USE_LEDA
#if defined(LEDA_NAMESPACE)
 std::cout << "leda::print_statistics() is missing in the free edition" << std::endl;
//  leda::print_statistics();
#else
//  print_statistics();
#endif
#endif
  destroy(V,n);
  return 0;
}
