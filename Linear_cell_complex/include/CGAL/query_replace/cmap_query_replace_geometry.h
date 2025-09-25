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
#ifndef CMAP_QUERY_REPLACE_GEOMETRY_H
#define CMAP_QUERY_REPLACE_GEOMETRY_H

#include <unordered_map>
#include <utility>
#include <vector>
#include <tuple>
#include <queue>

#include <CGAL/Kernel_traits.h>

#include "cmap_3close_cc.h"
#include "cmap_signature.h"
#include "lcc_geometry_transformation.h"
#include "lcc_read_depending_extension.h"

template<class LCC>
class Pattern_substituer;

template<class LCC, unsigned int>
class Pattern;

template<class LCC>
using Dart_mapping=std::unordered_map<typename LCC::Dart_handle,
                                      typename LCC::Dart_handle>;

///////////////////////////////////////////////////////////////////////////////
template<class LCC, unsigned int type>
class Barycentric_coord
{};

template<class LCC>
class Barycentric_coord<LCC, 1>
{
public:
  using Dart_handle=typename LCC::Dart_handle;

  void display(LCC& lcc)
  {
    std::cout<<lcc.point(m_dart)<<": ";
    for(auto& it: m_coords)
    { std::cout<<"["<<"  "<<lcc.point(std::get<0>(it))<<": "<<std::get<1>(it)
               <<", "<<std::get<2>(it)<<", "<<std::get<3>(it)<<"] "; }
    std::cout<<std::endl;
  }

  Dart_handle m_dart; // dart of an inner vertex
  /// barycentric coords of this inner vertex for each border vertex
  std::vector<std::tuple<Dart_handle, double, double, double>> m_coords;
};
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
class Barycentric_coord<LCC, 3>
{
public:
  using Dart_handle=typename LCC::Dart_handle;

  void display(LCC& lcc)
  {
    std::cout<<lcc.point(m_dart)<<": ";
    for(auto& it: m_coords)
    { std::cout<<"["<<"  "<<lcc.point(std::get<0>(it))<<": "<<std::get<1>(it)<<", "
               <<", "<<std::get<2>(it)<<", "<<std::get<3>(it)<<", "
               <<std::get<4>(it)<<"] "; }
    std::cout<<std::endl;
  }

  Dart_handle m_dart;
  std::vector<std::tuple<Dart_handle, double, double, double, double>> m_coords;
};
///////////////////////////////////////////////////////////////////////////////
template<typename Point>
bool compute_alpha_beta_gamma_of_point(const Point& a, const Point& b,
                                       const Point& c, const Point& p,
                                       double& alpha, double& beta, double& gamma)
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
template<typename LCC>
typename LCC::Point
compute_point_2D(LCC& lcc,
                 Dart_mapping<LCC>& links_from_pattern_to_face,
                 Dart_mapping<LCC>& pattern_to_global,
                 Barycentric_coord<LCC, 1>& m_barycentric_coords)
{
  typename LCC::Point p;
  typename LCC::Vector res=CGAL::NULL_VECTOR;
  typename LCC::Dart_handle cur, dh1, dh2;
  std::size_t nb=0;
  typename LCC::Dart_handle firstdh=
      links_from_pattern_to_face[pattern_to_global
      [std::get<0>(m_barycentric_coords.m_coords.front())]];
  typename LCC::Point bary2=lcc.template barycenter<2>(firstdh);
  for(std::tuple<typename LCC::Dart_handle, double, double, double>& e:
      m_barycentric_coords.m_coords)
  {
    assert(pattern_to_global.find(std::get<0>(e))!=pattern_to_global.end());
    assert(links_from_pattern_to_face.find(pattern_to_global[std::get<0>(e)])
        !=links_from_pattern_to_face.end());
    cur=links_from_pattern_to_face[pattern_to_global[std::get<0>(e)]];
    dh1=lcc.other_extremity(cur);
    if(compute_point_from_alpha_beta_gamma(lcc.point(cur), lcc.point(dh1),
                                           bary2, std::get<1>(e),
                                           std::get<2>(e), std::get<3>(e), p))
    {
      res+=typename LCC::Vector(p.x(), p.y(), p.z());
      ++nb;
    }
  }
  assert(nb>0);
  return typename LCC::Point(res.x()/nb, res.y()/nb, res.z()/nb);
}
///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
typename LCC::Point
compute_point_3D(LCC& lcc,
                 Dart_mapping<LCC>& links_from_pattern_to_volume,
                 Dart_mapping<LCC>& pattern_to_global,
                 Barycentric_coord<LCC, 3>& m_barycentric_coords)
{
  typename LCC::Point p;
  typename LCC::Vector res=CGAL::NULL_VECTOR;
  typename LCC::Dart_handle cur, dh1, dh2, dh3;
  std::size_t nb=0;
  for(std::tuple<typename LCC::Dart_handle, double, double, double, double>& e:
      m_barycentric_coords.m_coords)
  {
    assert(pattern_to_global.find(std::get<0>(e))!=pattern_to_global.end());
    assert(links_from_pattern_to_volume.find(pattern_to_global[std::get<0>(e)])
        !=links_from_pattern_to_volume.end());
    cur=links_from_pattern_to_volume[pattern_to_global[std::get<0>(e)]];
    dh1=lcc.template beta<0>(cur);
    dh2=lcc.other_extremity(cur);
    dh3=lcc.template beta<2,1,2>(cur);
    if(compute_point_from_alpha_beta_gamma_delta(lcc.point(cur), lcc.point(dh1),
                                                 lcc.point(dh2), lcc.point(dh3),
                                                 std::get<1>(e), std::get<2>(e),
                                                 std::get<3>(e), std::get<4>(e),
                                                 p))
    {
      res+=typename LCC::Vector(p.x(), p.y(), p.z());
      ++nb;
    }
  }
  assert(nb>0);
  return typename LCC::Point(res.x()/nb, res.y()/nb, res.z()/nb);
}
///////////////////////////////////////////////////////////////////////////////
template<typename LCC>
typename LCC::Point
compute_point_3D_v2(LCC& lcc,
                    Dart_mapping<LCC>& links_from_pattern_to_volume,
                    Dart_mapping<LCC>& pattern_to_global,
                    Barycentric_coord<LCC, 3>& m_barycentric_coords)
{
  assert(!m_barycentric_coords.m_coords.empty());
  typename LCC::Point p;
  typename LCC::Vector res=CGAL::NULL_VECTOR;
  typename LCC::Dart_handle cur, dh1, dh2, dh3;
  std::size_t nb=0;
  typename LCC::Dart_handle firstdh=
      links_from_pattern_to_volume[pattern_to_global
      [std::get<0>(m_barycentric_coords.m_coords.front())]];
  typename LCC::Point bary3=lcc.template barycenter<3>(firstdh);
  for(std::tuple<typename LCC::Dart_handle, double, double, double, double>& e:
      m_barycentric_coords.m_coords)
  {
    assert(pattern_to_global.find(std::get<0>(e))!=pattern_to_global.end());
    assert(links_from_pattern_to_volume.find(pattern_to_global[std::get<0>(e)])
        !=links_from_pattern_to_volume.end());
    cur=links_from_pattern_to_volume[pattern_to_global[std::get<0>(e)]];
    dh1=lcc.other_extremity(cur);

    // TODO avoid to recompute barycenters several times (?)
    if(compute_point_from_alpha_beta_gamma_delta(lcc.point(cur), lcc.point(dh1),
                                                 lcc.template barycenter<2>(cur),
                                                 bary3,
                                                 std::get<1>(e), std::get<2>(e),
                                                 std::get<3>(e), std::get<4>(e),
                                                 p))
    {
      res+=typename LCC::Vector(p.x(), p.y(), p.z());
      ++nb;
    }
  }
  assert(nb>0);
  return typename LCC::Point(res.x()/nb, res.y()/nb, res.z()/nb);
}
///////////////////////////////////////////////////////////////////////////////
/// Transform the geometry of the fpattern according to the geometry of the
/// target.
template<typename LCC>
void transform_geometry_of_fpattern(LCC& lcc,
                                    Dart_mapping<LCC>& links_from_pattern_to_face,
                                    Dart_mapping<LCC>& pattern_to_global,
                                    Pattern<LCC, 1>& pattern)
{
  for(Barycentric_coord<LCC, 1>& inner: pattern.barycentric_coords())
  {
    assert(pattern_to_global.find(inner.m_dart)!=pattern_to_global.end());
    typename LCC::Dart_handle res=pattern_to_global[inner.m_dart];
    // TODO avoid to recompute barycenters several times (?)
    lcc.point(res)=compute_point_2D(lcc,
                                    links_from_pattern_to_face,
                                    pattern_to_global,
                                    inner);
  }
}
///////////////////////////////////////////////////////////////////////////////
/// Transform the geometry of the spattern according to the geometry of the
/// target. For now same method than transform_geometry_of_fpattern
template<typename LCC>
void transform_geometry_of_spattern(LCC& lcc,
                                    Dart_mapping<LCC>& links_from_pattern_to_face,
                                    Dart_mapping<LCC>& pattern_to_global,
                                    Pattern<LCC, 2>& pattern)
{
  for(Barycentric_coord<LCC, 1>& inner: pattern.barycentric_coords())
  {
    assert(pattern_to_global.find(inner.m_dart)!=pattern_to_global.end());
    typename LCC::Dart_handle res=pattern_to_global[inner.m_dart];
    // TODO avoid to recompute barycenters several times (?)
    lcc.point(res)=compute_point_2D(lcc,
                                    links_from_pattern_to_face,
                                    pattern_to_global,
                                    inner);
  }
}
////////////////////////////////////////////////////////////////////////////////
/// Transform the geometry of the vpattern according to the geometry of the
/// target. Mark the dart of the external faces of the pattern.
/// For now simple solution that does not work for any pattern. TODO better?
template<typename LCC>
void transform_geometry_of_vpattern(LCC& lcc,
                                    Dart_mapping<LCC>& links_from_pattern_to_volume,
                                    Dart_mapping<LCC>& pattern_to_global,
                                    Pattern<LCC, 3>& pattern)
{
  for(Barycentric_coord<LCC, 3>& inner: pattern.barycentric_coords())
  {
    assert(pattern_to_global.find(inner.m_dart)!=pattern_to_global.end());
    typename LCC::Dart_handle res=pattern_to_global[inner.m_dart];
    // std::cout<<"[transform_geometry_of_vpattern] "<<lcc.point()<<" -> before "
    //          <<lcc.point()<<" and  after ";
    /* lcc.point(res)=compute_point_3D(lcc,
                                    links_from_pattern_to_volume,
                                    pattern_to_global,
                                    inner); */
    lcc.point(res)=compute_point_3D_v2(lcc,
                                       links_from_pattern_to_volume,
                                       pattern_to_global,
                                       inner);
    // std::cout<<lcc.point()<<std::endl;
  }
}
///////////////////////////////////////////////////////////////////////////////
template<class LCC, unsigned int type> // type==1 for face, 2 for surface, 3 for volume
class Pattern;
////////////////////////////////////////////////////////////////////////////////
template<class LCC>
class Pattern<LCC, 1> // Face pattern
{
  friend class Pattern_substituer<LCC>;

  using Dart_handle=typename LCC::Dart_handle;
  using size_type=typename LCC::size_type;
  using Point=typename LCC::Point;
  using Vector=typename LCC::Vector;
public:
  Pattern(): m_mark_to_preserve(LCC::INVALID_MARK)
  {}

  LCC& lcc()
  { return m_lcc; }

  size_type reserve_mark_to_preserve()
  {
    if(m_mark_to_preserve==LCC::INVALID_MARK)
    { m_mark_to_preserve=m_lcc.get_new_mark(); }
    return m_mark_to_preserve;
  }

  size_type mark_to_preserve() const
  { return m_mark_to_preserve; }

  std::vector<Barycentric_coord<LCC, 1>>& barycentric_coords()
  { return m_barycentric_coords; }

  void compute_barycentric_coord()
  {
    auto mark_vertices=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end(); it!=itend; ++it)
    {
      if(!m_lcc.is_marked(it, mark_vertices))
      {
        if(m_lcc.template is_free<2>(it)) // not an inner vertex
        { m_lcc.template mark_cell<0>(it, mark_vertices); }
      }
    }

    typename LCC::Point bary2=CGAL::ORIGIN;
    std::vector<Dart_handle> boundary_darts;
    std::vector<Dart_handle> inner_vertices;
    std::size_t nb1=0;
    auto vertex_treated=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end(); it!=itend; ++it)
    {
      if(m_lcc.template is_free<2>(it))
      { boundary_darts.push_back(it); }
      if(!m_lcc.is_marked(it, vertex_treated))
      { // TODO we can improve the traversal of cells (regroup, use basic it...)
        m_lcc.template mark_cell<0>(it, vertex_treated);
        if(!m_lcc.is_marked(it, mark_vertices))
        { // Here it is incident to an inner vertex
          m_lcc.template mark_cell<0>(it, mark_vertices);
          inner_vertices.push_back(it);
        }
        else // Here is is incident to a vertex of the boundary
        {
          const Point& p=m_lcc.point(it);
          bary2=Point(bary2.x()+p.x(), bary2.y()+p.y(), bary2.z()+p.z());
          ++nb1;
        }
      }
    }
    m_lcc.free_mark(mark_vertices);
    m_lcc.free_mark(vertex_treated);
    assert(nb1>0);
    bary2=Point(bary2.x()/nb1, bary2.y()/nb1, bary2.z()/nb1);

    nb1=0;
    m_barycentric_coords.resize(inner_vertices.size());
    for(auto it: inner_vertices)
    {
      m_barycentric_coords[nb1].m_coords.reserve(boundary_darts.size());
      m_barycentric_coords[nb1].m_dart=it;
      ++nb1;
    }

    double alpha, beta, gamma;
    for(auto& itd: boundary_darts)
    {
      const Point& p0=m_lcc.point(itd);
      const Point& p1=m_lcc.point(m_lcc.other_extremity(itd));
      nb1=0;
      for(auto& it: inner_vertices)
      {
        if(compute_alpha_beta_gamma_of_point(p0, p1, bary2, m_lcc.point(it),
                                             alpha, beta, gamma))
        {
          m_barycentric_coords[nb1].m_coords.push_back
              (std::make_tuple(itd, alpha, beta, gamma));
        ++nb1;
        }
      }
    }
  }

  void display()
  {
    for(auto& it: m_barycentric_coords)
    { it.display(m_lcc); }
  }

protected:
  LCC m_lcc;
  typename LCC::size_type m_mark_to_preserve;
  /// For each inner point, its barycentric coordinates for each external point
  std::vector<Barycentric_coord<LCC, 1>> m_barycentric_coords;
};
///////////////////////////////////////////////////////////////////////////////
template<class LCC>
class Pattern<LCC, 2> // Surfacic pattern
{
  friend class Pattern_substituer<LCC>;

  using Dart_handle=typename LCC::Dart_handle;
  using size_type=typename LCC::size_type;
  using Point=typename LCC::Point;
  using Vector=typename LCC::Vector;
public:
  Pattern():
    m_mark_faceborder(m_lcc.get_new_mark()),
    m_mark_to_preserve(LCC::INVALID_MARK)
  {}

  LCC& lcc()
  { return m_lcc; }

  size_type reserve_mark_to_preserve()
  {
    if(m_mark_to_preserve==LCC::INVALID_MARK)
    { m_mark_to_preserve=m_lcc.get_new_mark(); }
    return m_mark_to_preserve;
  }

  size_type mark_to_preserve() const
  { return m_mark_to_preserve; }

  // same barycentric coords than for face => 1
  std::vector<Barycentric_coord<LCC, 1>>& barycentric_coords()
  { return m_barycentric_coords; }

  void compute_barycentric_coord()
  {
    // Mark vertices that belong to a face border.
    auto border_vertex=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end(); it!=itend; ++it)
    {
      if(!m_lcc.is_marked(it, border_vertex))
      {
        if(m_lcc.is_marked(it, m_mark_faceborder)) // not an inner vertex
        { m_lcc.template mark_cell<0>(it, border_vertex); }
      }
    }

    typename LCC::Point bary2=CGAL::ORIGIN;
    std::vector<Dart_handle> boundary_darts;
    std::vector<Dart_handle> inner_vertices;
    std::size_t nb1=0, new_index=0;
    Dart_handle cur, other;
    auto treated=m_lcc.get_new_mark();
    auto vertex_treated=m_lcc.get_new_mark();
    std::queue<Dart_handle> to_treat;
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end(); it!=itend; ++it)
    {
      if(!m_lcc.is_marked(it, treated))
      {
        boundary_darts.clear();
        inner_vertices.clear();
        nb1=0;
        bary2=CGAL::ORIGIN;

        // Here we iterate through the cc of faces inside a same cycle of edges
        // marked by m_mark_faceborder
        to_treat.push(it);
        m_lcc.mark(it, treated);
        while(!to_treat.empty())
        {
          cur=to_treat.front();
          to_treat.pop();

          if(m_lcc.is_marked(cur, m_mark_faceborder))
          {  boundary_darts.push_back(cur); }

          if(!m_lcc.is_marked(cur, vertex_treated))
          { // TODO we can improve the traversal of cells (regroup, use basic it...)
            m_lcc.template mark_cell<0>(cur, vertex_treated);
            if(!m_lcc.is_marked(cur, border_vertex))
            { // Here it is incident to an inner vertex
              inner_vertices.push_back(cur);
            }
            else // Here is is incident to a vertex of the boundary
            {
              const Point& p=m_lcc.point(cur);
              bary2=Point(bary2.x()+p.x(), bary2.y()+p.y(), bary2.z()+p.z());
              ++nb1;
            }
          }

          other=m_lcc.template beta<1>(cur);
          if(!m_lcc.is_marked(other, treated))
          {
            to_treat.push(other);
            m_lcc.mark(other, treated);
          }

          if(!m_lcc.is_marked(cur, m_mark_faceborder))
          {
            other=m_lcc.template beta<2>(cur);
            if(!m_lcc.is_marked(other, treated))
            {
              to_treat.push(other);
              m_lcc.mark(other, treated);
            }
          }
        }

        assert(nb1>0);
        // Now compute the barycentric coordinates of inner vertices
        if(inner_vertices.size()>0)
        {
          bary2=Point(bary2.x()/nb1, bary2.y()/nb1, bary2.z()/nb1);

          new_index=m_barycentric_coords.size();
          nb1=new_index;
          m_barycentric_coords.resize(m_barycentric_coords.size()+
                                      inner_vertices.size());
          for(auto iti: inner_vertices)
          {
            m_barycentric_coords[nb1].m_coords.reserve(boundary_darts.size());
            m_barycentric_coords[nb1].m_dart=iti;
            ++nb1;
            if(m_lcc.is_marked(iti, vertex_treated))
            { m_lcc.template unmark_cell<0>(iti, vertex_treated);}
          }

          double alpha, beta, gamma;
          for(auto& itd: boundary_darts)
          {
            const Point& p0=m_lcc.point(itd);
            const Point& p1=m_lcc.point(m_lcc.other_extremity(itd));
            nb1=new_index;
            for(auto& iti: inner_vertices)
            {
              if(compute_alpha_beta_gamma_of_point(p0, p1, bary2, m_lcc.point(iti),
                                                   alpha, beta, gamma))
              {
                m_barycentric_coords[nb1].m_coords.push_back
                    (std::make_tuple(itd, alpha, beta, gamma));
                ++nb1;
              }
            }
            if(m_lcc.is_marked(itd, vertex_treated))
            { m_lcc.template unmark_cell<0>(itd, vertex_treated);}
          }
        }
        else
        {
          for(auto& itd: boundary_darts)
          {
            if(m_lcc.is_marked(itd, vertex_treated))
            { m_lcc.template unmark_cell<0>(itd, vertex_treated);}
          }
        }
        assert(m_lcc.is_whole_map_unmarked(vertex_treated));
      }
    }
    assert(m_lcc.is_whole_map_marked(treated));
    m_lcc.free_mark(treated);
    m_lcc.free_mark(border_vertex);
    m_lcc.free_mark(vertex_treated);
 }

  void display()
  {
    for(auto& it: m_barycentric_coords)
    { it.display(m_lcc); }
  }

protected:
  LCC m_lcc;
  typename LCC::size_type m_mark_faceborder;
  typename LCC::size_type m_mark_to_preserve;
  /// For each inner point, its barycentric coordinates for each external point
  std::vector<Barycentric_coord<LCC, 1>> m_barycentric_coords;
};
////////////////////////////////////////////////////////////////////////////////
template<class LCC>
class Pattern<LCC, 3> // Volumic pattern
{
  friend class Pattern_substituer<LCC>;

  using Dart_handle=typename LCC::Dart_handle;
  using size_type=typename LCC::size_type;
  using Point=typename LCC::Point;
  using Vector=typename LCC::Vector;
public:
  Pattern(): m_mark_to_preserve(LCC::INVALID_MARK)
  {}

  LCC& lcc()
  { return m_lcc; }

  size_type reserve_mark_to_preserve()
  {
    if(m_mark_to_preserve==LCC::INVALID_MARK)
    { m_mark_to_preserve=m_lcc.get_new_mark(); }
    return m_mark_to_preserve;
  }

  size_type mark_to_preserve() const
  { return m_mark_to_preserve; }

  std::vector<Barycentric_coord<LCC, 3>>& barycentric_coords()
  { return m_barycentric_coords; }

  void compute_barycentric_coord()
  {
    auto mark_vertices=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end(); it!=itend; ++it)
    {
      if(!m_lcc.is_marked(it, mark_vertices))
      {
        if(m_lcc.template is_free<3>(it)) // not an inner vertex
        { m_lcc.template mark_cell<0>(it, mark_vertices); }
      }
    }

    typename LCC::Point bary3=CGAL::ORIGIN;
    std::vector<Dart_handle> boundary_darts;
    std::vector<Dart_handle> inner_vertices;
    std::size_t nb1=0;
    auto vertex_treated=m_lcc.get_new_mark();
    for(auto it=m_lcc.darts().begin(), itend=m_lcc.darts().end(); it!=itend; ++it)
    {
      if(m_lcc.template is_free<3>(it))
      { boundary_darts.push_back(it); }
      if(!m_lcc.is_marked(it, vertex_treated))
      { // TODO we can improve the traversal of cells (regroup, use basic it...)
        m_lcc.template mark_cell<0>(it, vertex_treated);
        if(!m_lcc.is_marked(it, mark_vertices))
        { // Here it is incident to an inner vertex
          m_lcc.template mark_cell<0>(it, mark_vertices);
          inner_vertices.push_back(it);
        }
        else // Here is is incident to a vertex of the boundary
        {
          const Point& p=m_lcc.point(it);
          bary3=Point(bary3.x()+p.x(), bary3.y()+p.y(), bary3.z()+p.z());
          ++nb1;
        }
      }
    }
    m_lcc.free_mark(mark_vertices);
    m_lcc.free_mark(vertex_treated);
    assert(nb1>0);
    bary3=Point(bary3.x()/nb1, bary3.y()/nb1, bary3.z()/nb1);

    nb1=0;
    m_barycentric_coords.resize(inner_vertices.size());
    for(auto it: inner_vertices)
    {
      m_barycentric_coords[nb1].m_coords.reserve(boundary_darts.size());
      m_barycentric_coords[nb1].m_dart=it;
      ++nb1;
    }

    double alpha, beta, gamma, delta;
    for(auto& itd: boundary_darts)
    {
      const Point& p0=m_lcc.point(itd);
      const Point& p1=m_lcc.point(m_lcc.other_extremity(itd));
      const Point p2=m_lcc.template barycenter<2>(itd);
      nb1=0;
      for(auto& it: inner_vertices)
      { // TODO we can test if tetra (p0, p1, p2, bary3) before to enter in the loop
        if(compute_alpha_beta_gamma_delta_of_point(p0, p1, p2, bary3,
                                                   m_lcc.point(it),
                                                   alpha, beta, gamma, delta))
        {
          m_barycentric_coords[nb1].m_coords.push_back
              (std::make_tuple(itd, alpha, beta, gamma, delta));
          ++nb1;
        }
      }
    }
  }

  void display()
  {
    for(auto& it: m_barycentric_coords)
    { it.display(m_lcc); }
  }

protected:
  LCC m_lcc;
  size_type m_mark_to_preserve;
  /// For each inner point, its barycentric coordinates for each external point
  std::vector<Barycentric_coord<LCC, 3>> m_barycentric_coords;
};
////////////////////////////////////////////////////////////////////////////////
#endif // CMAP_QUERY_REPLACE_GEOMETRY_H
