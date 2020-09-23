// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
// $Date$

// Author(s)     :  Manuel Caroli <Manuel.Caroli@sophia.inria.fr>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Simple_cartesian.h>

#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>

#include <CGAL/Random.h>

typedef CGAL::Simple_cartesian<CGAL::Quotient<CGAL::MP_Float> > EK;
typedef CGAL::Exact_predicates_inexact_constructions_kernel     FK;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<EK>             ETraits;
typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<FK>             FTraits;

typedef ETraits::Point_3              EPoint;
typedef ETraits::Iso_cuboid_3         EIso_cuboid;
typedef ETraits::Periodic_3_offset_3  EOffset;
typedef FTraits::Point_3              FPoint;
typedef FTraits::Iso_cuboid_3         FIso_cuboid;
typedef FTraits::Periodic_3_offset_3  FOffset;

typedef ETraits::Orientation_3 EOrientation_3;
typedef FTraits::Orientation_3 FOrientation_3;
typedef ETraits::Side_of_oriented_sphere_3 ESide_of_oriented_sphere_3;
typedef FTraits::Side_of_oriented_sphere_3 FSide_of_oriented_sphere_3;

typedef FPoint  Point;
typedef FOffset Offset;

CGAL::Random *r;

double rand_base()
{
  return r->get_double(0, 1);
}

// Random double almost in [0;1].
double my_rand_dbl()
{
  // Ensure 53 random bits, not 48.
  return rand_base() + rand_base()/1024;
}

// Random point in unit cube.
Point my_rand_p3()
{
  double x = my_rand_dbl(), y = my_rand_dbl(), z = my_rand_dbl();
  return Point(x, y, z);
}

// Random int in [0;256).
int my_rand_int(int imin, int imax)
{
  return r->get_int(imin, imax+1);
}

// Random offset
Offset my_rand_o3(int imin, int imax)
{
  int x = my_rand_int(imin,imax);
  int y = my_rand_int(imin,imax);
  int z = my_rand_int(imin,imax);
  return Offset(x, y, z);
}

// Random point-offset pair on sphere of radius 1 and center (1,1,1).
std::pair<Point,Offset> sphere_rand_pp3()
{
  double dx = my_rand_dbl();
  double dy = my_rand_dbl();
  double dz = my_rand_dbl();
  double n = std::sqrt(dx*dx + dy*dy + dz*dz);
  double x = dx/n, y = dy/n, z = dz/n;

  int ix = my_rand_int(0,1);
  int iy = my_rand_int(0,1);
  int iz = my_rand_int(0,1);
  Offset o(ix,iy,iz);

  Point p(((ix==0)?1-x:x),((iy==0)?1-y:y),((iz==0)?1-z:z));
  return std::make_pair(p,o);
}

EPoint ept(Point &p) {
  return EPoint(p.x(),p.y(),p.z());
}

EOffset eoff(Offset &o) {
  return EOffset(o.x(),o.y(),o.z());
}

// Perturbation with given maximum relative epsilon.
void perturb(Point &p, double rel_eps)
{
  double x = p.x()*(1+rand_base()*rel_eps);
  double y = p.y()*(1+rand_base()*rel_eps);
  double z = p.z()*(1+rand_base()*rel_eps);
  p = Point(x, y, z);
}

// Pick a random point on the triangle [p, q, r].
std::pair<Point, Offset> pick_coplanar(
    const Point &p, const Point &q, const Point &r,
    const Offset &po, const Offset &qo, const Offset& ro)
{
  // s = p + (p-q)*my_rand() + (p-r)*my_rand();  (almost)

  double r1 = my_rand_dbl(), r2 = my_rand_dbl();
  double x = p.x()+po.x() + (q.x()+qo.x() - (p.x()+po.x()))*r1
                          + (r.x()+ro.x() - (p.x()+po.x()))*r2;
  double y = p.y()+po.y() + (q.y()+qo.y() - (p.y()+po.y()))*r1
                          + (r.y()+ro.y() - (p.y()+po.y()))*r2;
  double z = p.z()+po.z() + (q.z()+qo.z() - (p.z()+po.z()))*r1
                          + (r.z()+ro.z() - (p.z()+po.z()))*r2;

  int ix = static_cast<int>(floor(x));
  int iy = static_cast<int>(floor(y));
  int iz = static_cast<int>(floor(z));

  return std::make_pair(Point(x-floor(x),y-floor(y),z-floor(z)),
      Offset(ix,iy,iz));
}

void
test_orientation_3(EOrientation_3 & eorient, FOrientation_3 & forient)
{
  // First test with random points.
  Point fp = my_rand_p3();
  Point fq = my_rand_p3();
  Point fr = my_rand_p3();
  Point fs = my_rand_p3();
  EPoint ep = ept(fp);
  EPoint eq = ept(fq);
  EPoint er = ept(fr);
  EPoint es = ept(fs);

  Offset fpo = my_rand_o3(0,255);
  Offset fqo = my_rand_o3(0,255);
  Offset fro = my_rand_o3(0,255);
  Offset fso = my_rand_o3(0,255);
  Offset epo = eoff(fpo);
  Offset eqo = eoff(fqo);
  Offset ero = eoff(fro);
  Offset eso = eoff(fso);

  assert(forient(fp,fq,fr,fs,fpo,fqo,fro,fso)
      == eorient(ep,eq,er,es,epo,eqo,ero,eso));
  assert(forient(fq,fr,fs,fp,fqo,fro,fso,fpo)
      == eorient(eq,er,es,ep,eqo,ero,eso,epo));
  assert(forient(fr,fs,fp,fq,fro,fso,fpo,fqo)
      == eorient(er,es,ep,eq,ero,eso,epo,eqo));
  assert(forient(fp,fq,fr,fs,fpo,fqo,fro,fso)
      == eorient(ep,eq,er,es,epo,eqo,ero,eso));

  // Then with coplanar points (up to roundoff errors).
  //s = p + (p-q)*my_rand() + (p-r)*my_rand();
  std::pair<Point, Offset> fspo = pick_coplanar(fp, fq, fr, fpo, fqo, fro);
  fs = fspo.first;
  es = ept(fs);
  fso = fspo.second;
  eso = eoff(fso);

  assert(forient(fp,fq,fr,fs,fpo,fqo,fro,fso)
      == eorient(ep,eq,er,es,epo,eqo,ero,eso));
  assert(forient(fq,fr,fs,fp,fqo,fro,fso,fpo)
      == eorient(eq,er,es,ep,eqo,ero,eso,epo));
  assert(forient(fr,fs,fp,fq,fro,fso,fpo,fqo)
      == eorient(er,es,ep,eq,ero,eso,epo,eqo));
  assert(forient(fp,fq,fr,fs,fpo,fqo,fro,fso)
      == eorient(ep,eq,er,es,epo,eqo,ero,eso));

  // Then with some perturbation.
  perturb(fs, 1.0/(1<<20)/(1<<20)); // 2^-40
  es = ept(fs);

  assert(forient(fp,fq,fr,fs,fpo,fqo,fro,fso)
      == eorient(ep,eq,er,es,epo,eqo,ero,eso));
  assert(forient(fq,fr,fs,fp,fqo,fro,fso,fpo)
      == eorient(eq,er,es,ep,eqo,ero,eso,epo));
  assert(forient(fr,fs,fp,fq,fro,fso,fpo,fqo)
      == eorient(er,es,ep,eq,ero,eso,epo,eqo));
  assert(forient(fp,fq,fr,fs,fpo,fqo,fro,fso)
      == eorient(ep,eq,er,es,epo,eqo,ero,eso));
}

void
test_side_of_oriented_sphere_3(ESide_of_oriented_sphere_3 & esoos,
    FSide_of_oriented_sphere_3 & fsoos)
{
  // First test with random points.
  Point fp = my_rand_p3();
  Point fq = my_rand_p3();
  Point fr = my_rand_p3();
  Point fs = my_rand_p3();
  Point ft = my_rand_p3();
  EPoint ep = ept(fp);
  EPoint eq = ept(fq);
  EPoint er = ept(fr);
  EPoint es = ept(fs);
  EPoint et = ept(ft);

  Offset fpo = my_rand_o3(0,255);
  Offset fqo = my_rand_o3(0,255);
  Offset fro = my_rand_o3(0,255);
  Offset fso = my_rand_o3(0,255);
  Offset fto = my_rand_o3(0,255);
  Offset epo = eoff(fpo);
  Offset eqo = eoff(fqo);
  Offset ero = eoff(fro);
  Offset eso = eoff(fso);
  Offset eto = eoff(fto);

  assert(fsoos(fp,fq,fr,fs,ft,fpo,fqo,fro,fso,fto)
      == esoos(ep,eq,er,es,et,epo,eqo,ero,eso,eto));
  assert(fsoos(fq,fr,fs,ft,fp,fqo,fro,fso,fto,fpo)
      == esoos(eq,er,es,et,ep,eqo,ero,eso,eto,epo));
  assert(fsoos(fr,fs,ft,fp,fq,fro,fso,fto,fpo,fqo)
      == esoos(er,es,et,ep,eq,ero,eso,eto,epo,eqo));
  assert(fsoos(fs,ft,fp,fq,fr,fso,fto,fpo,fqo,fro)
      == esoos(es,et,ep,eq,er,eso,eto,epo,eqo,ero));
  assert(fsoos(ft,fs,fp,fq,fr,fto,fso,fpo,fqo,fro)
      == esoos(et,es,ep,eq,er,eto,eso,epo,eqo,ero));

  // Then with cospherical points (up to roundoff errors).
  std::pair<Point,Offset> fppo = sphere_rand_pp3();
  std::pair<Point,Offset> fqpo = sphere_rand_pp3();
  std::pair<Point,Offset> frpo = sphere_rand_pp3();
  std::pair<Point,Offset> fspo = sphere_rand_pp3();
  std::pair<Point,Offset> ftpo = sphere_rand_pp3();
  fp = fppo.first;
  fq = fqpo.first;
  fr = frpo.first;
  fs = fspo.first;
  ft = ftpo.first;
  ep = ept(fp);
  eq = ept(fq);
  er = ept(fr);
  es = ept(fs);
  et = ept(ft);
  fpo = fppo.second;
  fqo = fqpo.second;
  fro = frpo.second;
  fso = fspo.second;
  fto = ftpo.second;
  epo = eoff(fpo);
  eqo = eoff(fqo);
  ero = eoff(fro);
  eso = eoff(fso);
  eto = eoff(fto);

  assert(fsoos(fp,fq,fr,fs,ft,fpo,fqo,fro,fso,fto)
      == esoos(ep,eq,er,es,et,epo,eqo,ero,eso,eto));
  assert(fsoos(fq,fr,fs,ft,fp,fqo,fro,fso,fto,fpo)
      == esoos(eq,er,es,et,ep,eqo,ero,eso,eto,epo));
  assert(fsoos(fr,fs,ft,fp,fq,fro,fso,fto,fpo,fqo)
      == esoos(er,es,et,ep,eq,ero,eso,eto,epo,eqo));
  assert(fsoos(fs,ft,fp,fq,fr,fso,fto,fpo,fqo,fro)
      == esoos(es,et,ep,eq,er,eso,eto,epo,eqo,ero));
  assert(fsoos(ft,fs,fp,fq,fr,fto,fso,fpo,fqo,fro)
      == esoos(et,es,ep,eq,er,eto,eso,epo,eqo,ero));

  // Then with some perturbation.
  perturb(fr, 1.0/(1<<25)/(1<<20)); // 2^-45
  er = ept(fr);

  assert(fsoos(fp,fq,fr,fs,ft,fpo,fqo,fro,fso,fto)
      == esoos(ep,eq,er,es,et,epo,eqo,ero,eso,eto));
  assert(fsoos(fq,fr,fs,ft,fp,fqo,fro,fso,fto,fpo)
      == esoos(eq,er,es,et,ep,eqo,ero,eso,eto,epo));
  assert(fsoos(fr,fs,ft,fp,fq,fro,fso,fto,fpo,fqo)
      == esoos(er,es,et,ep,eq,ero,eso,eto,epo,eqo));
  assert(fsoos(fs,ft,fp,fq,fr,fso,fto,fpo,fqo,fro)
      == esoos(es,et,ep,eq,er,eso,eto,epo,eqo,ero));
  assert(fsoos(ft,fs,fp,fq,fr,fto,fso,fpo,fqo,fro)
      == esoos(et,es,ep,eq,er,eto,eso,epo,eqo,ero));
}

void compute_epsilons()
{
  FOrientation_3::compute_epsilon();
  FSide_of_oriented_sphere_3::compute_epsilon();
}

int _test_periodic_3_static_filters()
{
  int loops = 2000;
  int seed  = CGAL::get_default_random().get_int(0, 1<<30);

  std::cout << "Initializing random generator with seed = " << seed
            << std::endl
            << "#loops = " << loops
            << std::endl;

  CGAL::Random rnd(seed);
  r = &rnd;

  std::cout.precision(20);
  std::cerr.precision(20);

  compute_epsilons();

  ETraits et;
  FTraits ft;
  et.set_domain(EIso_cuboid(0,0,0,1,1,1));
  ft.set_domain(FIso_cuboid(0,0,0,1,1,1));

  EOrientation_3 eorient = et.orientation_3_object();
  FOrientation_3 forient = ft.orientation_3_object();
  ESide_of_oriented_sphere_3 esoos = et.side_of_oriented_sphere_3_object();
  FSide_of_oriented_sphere_3 fsoos = ft.side_of_oriented_sphere_3_object();

  std::cout << "Testing statically filtered Orientation_3" << std::endl;
  for(int i=0; i<loops; ++i)
    test_orientation_3(eorient,forient);

  std::cout << "Testing statically filtered Side_of_oriented_sphere_3"
            << std::endl;
  for(int i=0; i<loops; ++i)
    test_side_of_oriented_sphere_3(esoos,fsoos);

  return 0;
}

