// Copyright (c) 2025 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
////////////////////////////////////////////////////////////////////////////////
#ifndef CGAL_LCC_BARYCENTRIC_COORD_H
#define CGAL_LCC_BARYCENTRIC_COORD_H

#include <iostream>
#include <tuple>
#include <vector>

#include <CGAL/Kernel_traits.h>
#include <CGAL/Kernel/global_functions.h>

namespace CGAL::internal
{
///////////////////////////////////////////////////////////////////////////////
template<class LCC, unsigned int type>
class Barycentric_coord
{};
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
class Barycentric_coord<LCC, 1>
{
public:
  using Dart_descriptor=typename LCC::Dart_descriptor;

  void display(LCC& lcc)
  {
    std::cout<<lcc.point(m_dart)<<": ";
    for(auto& it: m_coords)
    { std::cout<<"["<<"  "<<lcc.point(std::get<0>(it))<<": "<<std::get<1>(it)
               <<", "<<std::get<2>(it)<<", "<<std::get<3>(it)<<"] "; }
    std::cout<<std::endl;
  }

  Dart_descriptor m_dart; // dart of an inner vertex
  /// barycentric coords of this inner vertex for each border vertex
  std::vector<std::tuple<Dart_descriptor, double, double, double>> m_coords;
};
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
class Barycentric_coord<LCC, 3>
{
public:
  using Dart_descriptor=typename LCC::Dart_descriptor;

  void display(LCC& lcc)
  {
    std::cout<<lcc.point(m_dart)<<": ";
    for(auto& it: m_coords)
    { std::cout<<"["<<"  "<<lcc.point(std::get<0>(it))<<": "<<std::get<1>(it)
               <<", "<<std::get<2>(it)<<", "<<std::get<3>(it)<<", "
               <<std::get<4>(it)<<"] "; }
    std::cout<<std::endl;
  }

  Dart_descriptor m_dart;
  std::vector<std::tuple<Dart_descriptor, double, double, double, double>>
      m_coords;
};
///////////////////////////////////////////////////////////////////////////////
template<typename Point>
bool compute_alpha_beta_gamma_of_point(const Point& a, const Point& b,
                                       const Point& c, const Point& p,
                                       double& alpha, double& beta,
                                       double& gamma)
{
  typename CGAL::Kernel_traits<Point>::Kernel::Triangle_3 t(a, b, c);
  if(t.is_degenerate()) { return false; }

  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vap(a, p);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vbp(b, p);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vcp(c, p);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vab(a, b);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vac(a, c);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vca(c, a);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vbc(b, c);

  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3
      n=CGAL::cross_product(vab, vac);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3
      na=CGAL::cross_product(vbc, vbp);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3
      nb=CGAL::cross_product(vca, vcp);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3
      nc=CGAL::cross_product(vab, vap);

  alpha=(n*na)/(n*n);
  beta=(n*nb)/(n*n);
  gamma=(n*nc)/(n*n);
  return true;

 /* typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 v0(p0, p1);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 v1(p0, p2);

  double div=v0.y()*v1.z()-v1.y()*v0.z();
  if(div!=0)
  { beta=(v0.y()*(p.z()-p0.z())-(p.y()-p0.y())*v0.z())/(div); }
  else
  {
    div=v0.x()*v1.z()-v1.x()*v0.z();
    if(div!=0)
    { beta=(v0.x()*(p.z()-p0.z())-(p.x()-p0.x())*v0.z())/(div); }
    else
    {
      div=v0.y()*v1.x()-v1.y()*v0.x();
      assert(div!=0);
      beta=(v0.y()*(p.x()-p0.x())-(p.y()-p0.y())*v0.x())/(div);
    }
  }

  if(v0.x()!=0)
  { alpha=((p.x()-p0.x())-beta*v1.x())/v0.x(); }
  else if (v0.y()!=0)
  { alpha=((p.y()-p0.y())-beta*v1.y())/v0.y(); }
  else
  {
    assert(v1.z()!=0);
    alpha=((p.z()-p0.z())-beta*v1.z())/v0.z();
  }*/
}
///////////////////////////////////////////////////////////////////////////////
template<typename Point>
bool compute_point_from_alpha_beta_gamma(const Point& p0, const Point& p1,
                                         const Point& p2, double alpha,
                                         double beta, double gamma,
                                         Point& p)
{
  typename CGAL::Kernel_traits<Point>::Kernel::Triangle_3 t(p0, p1, p2);
  if(t.is_degenerate()) { return false; }
  p=Point(alpha*p0.x()+beta*p1.x()+gamma*p2.x(),
          alpha*p0.y()+beta*p1.y()+gamma*p2.y(),
          alpha*p0.z()+beta*p1.z()+gamma*p2.z());
  return true;

  /* typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 v0(p0, p1);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 v1(p0, p2);
  p=Point(p0.x()+alpha*v0.x()+beta*v1.x(),
          p0.y()+alpha*v0.y()+beta*v1.y(),
          p0.z()+alpha*v0.z()+beta*v1.z()); */
}
///////////////////////////////////////////////////////////////////////////////
template<typename Point>
bool compute_alpha_beta_gamma_delta_of_point(const Point& a, const Point& b,
                                             const Point& c, const Point& d,
                                             const Point& p, double& alpha,
                                             double& beta, double& gamma,
                                             double& delta)
{
  typename CGAL::Kernel_traits<Point>::Kernel::Tetrahedron_3 t(a, b, c, d);
  if(t.is_degenerate()) { return false; }

  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vap(a, p);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vbp(b, p);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vcp(c, p);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vdp(d, p);

  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vab(a, b);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vac(a, c);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vad(a, d);

  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vbc(b, c);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 vbd(b, d);

  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3
      temp=CGAL::cross_product(vac, vad);

  double va=(vbp*CGAL::cross_product(vbd, vbc));
  double vb=(vap*temp);
  double vc=(vap*CGAL::cross_product(vad, vab));
  double vd=(vap*CGAL::cross_product(vab, vac));
  double v=/* std::abs */((vab*temp));

  alpha=va/v;
  beta=vb/v;
  gamma=vc/v;
  delta=vd/v;
  return true;

  /* std::cout<<"[compute alpha...] "<<a<<" "<<b<<" "<<c<<" "<<d<<" "
           <<alpha<<" "<<beta<<" "<<gamma<<" "<<delta<<" -> "<<p<<std::endl; */

  /*
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 u(p0, p1);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 v(p0, p2);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 w(p0, p3);

  if((v.x()*u.y()-v.y()*u.x())!=0)
  {
    double part1=((p.x()*u.y()-p0.x()*u.y()-u.x()*p.y()+p0.y()*u.x())*
                  (v.z()*u.x()-v.x()*u.z()))/(v.x()*u.y()-v.y()*u.x());
    double part2=((w.y()*u.x()-w.x()*u.y())*(v.z()*u.x()-v.x()*u.z())+
                  (v.x()*u.y()-v.y()*u.x())*(w.z()*u.x()-w.x()*u.z()))/
                 (v.x()*u.y()-v.y()*u.x());
    assert(part2!=0);
    gamma=((p.z()*u.x()-p0.z()*u.x()-p.x()*u.z()+p0.x()*u.z())-part1)/part2;
  }
  else if((v.x()*u.z()-v.z()*u.x())!=0)
  {
    double part1=((p.x()*u.z()-p0.x()*u.z()-u.x()*p.z()+p0.z()*u.x())*
                  (v.y()*u.x()-v.x()*u.y()))/(v.x()*u.z()-v.z()*u.x());
    double part2=((w.z()*u.x()-w.x()*u.z())*(v.y()*u.x()-v.x()*u.y())+
                  (v.x()*u.z()-v.z()*u.x())*(w.y()*u.x()-w.x()*u.y()))/
                 (v.x()*u.z()-v.z()*u.x());
    assert(part2!=0);
    gamma=((p.y()*u.x()-p0.y()*u.x()-p.x()*u.y()+p0.x()*u.y())-part1)/part2;
  }
  else
  {
    assert((v.y()*u.z()-v.z()*u.y())!=0);
    double part1=((p.y()*u.z()-p0.y()*u.z()-u.y()*p.z()+p0.z()*u.y())*
                  (v.x()*u.y()-v.y()*u.x()))/(v.y()*u.z()-v.z()*u.y());
    double part2=((w.z()*u.y()-w.y()*u.z())*(v.x()*u.y()-v.y()*u.x())+
                  (v.y()*u.z()-v.z()*u.y())*(w.x()*u.y()-w.y()*u.x()))/
                 (v.y()*u.z()-v.z()*u.y());
    assert(part2!=0);
    gamma=((p.x()*u.y()-p0.x()*u.y()-p.y()*u.x()+p0.y()*u.x())-part1)/part2;
  }

  if((v.x()*u.y()-v.y()*u.x())!=0)
  {
    beta=((p.x()*u.y()-p0.x()*u.y()-u.x()*p.y()+p0.y()*u.x())/
          (v.x()*u.y()-v.y()*u.x()))+
         ((w.y()*u.x()-w.x()*u.y())/(v.x()*u.y()-v.y()*u.x()))*gamma;
  }
  else if((v.x()*u.z()-v.z()*u.x())!=0)
  {
    beta=((p.x()*u.z()-p0.x()*u.z()-u.x()*p.z()+p0.z()*u.x())/
          (v.x()*u.z()-v.z()*u.x()))+
         ((w.z()*u.x()-w.x()*u.z())/(v.x()*u.z()-v.z()*u.x()))*gamma;
  }
  else
  {
    assert((v.y()*u.z()-v.z()*u.y())!=0);
    beta=((p.y()*u.z()-p0.y()*u.z()-u.y()*p.z()+p0.z()*u.y())/
          (v.y()*u.z()-v.z()*u.y()))+
         ((w.z()*u.y()-w.y()*u.z())/(v.y()*u.z()-v.z()*u.y()))*gamma;
  }

  if(u.x()!=0)
  { alpha=(p.x()-p0.x()-beta*v.x()-gamma*w.x())/u.x(); }
  else if(u.y()!=0)
  { alpha=(p.y()-p0.y()-beta*v.y()-gamma*w.y())/u.y(); }
  else
  {
    assert(u.z()!=0);
    alpha=(p.z()-p0.z()-beta*v.z()-gamma*w.z())/u.z();
  }*/

  // TODO assertion true if we use epsilon comparison
  // assert(p.x()==alpha*a.x()+beta*b.x()+gamma*c.x()+delta*d.x());
  // assert(p.y()==alpha*a.y()+beta*b.y()+gamma*c.y()+delta*d.y());
  // assert(p.z()==alpha*a.z()+beta*b.z()+gamma*c.z()+delta*d.z());
}
///////////////////////////////////////////////////////////////////////////////
template<typename Point>
bool compute_point_from_alpha_beta_gamma_delta(const Point& p0, const Point& p1,
                                               const Point& p2, const Point& p3,
                                               double alpha, double beta,
                                               double gamma, double delta,
                                               Point& p)
{
  /* typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 u(p0, p1);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 v(p0, p2);
  typename CGAL::Kernel_traits<Point>::Kernel::Vector_3 w(p0, p3);
  p=Point(p0.x()+alpha*u.x()+beta*v.x()+gamma*w.x(),
          p0.y()+alpha*u.y()+beta*v.y()+gamma*w.y(),
          p0.z()+alpha*u.z()+beta*v.z()+gamma*w.z()); */
  typename CGAL::Kernel_traits<Point>::Kernel::Tetrahedron_3 t(p0, p1, p2, p3);
  if(t.is_degenerate()) { return false; }
  p=Point(alpha*p0.x()+beta*p1.x()+gamma*p2.x()+delta*p3.x(),
          alpha*p0.y()+beta*p1.y()+gamma*p2.y()+delta*p3.y(),
          alpha*p0.z()+beta*p1.z()+gamma*p2.z()+delta*p3.z());
  return true;
  /* std::cout<<"[compute point] "<<p0<<" "<<p1<<" "<<p2<<" "<<p3<<" "
           <<alpha<<" "<<beta<<" "<<gamma<<" "<<delta<<" -> "<<p<<std::endl; */
}
///////////////////////////////////////////////////////////////////////////////
}
#endif // LCC_BARYCENTRIC_COORD_H
