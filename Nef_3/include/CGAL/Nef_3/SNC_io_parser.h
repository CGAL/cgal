// Copyright (c) 1997-2002  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Michael Seel    <seel@mpi-sb.mpg.de>
//                 Miguel Granados <granados@mpi-sb.mpg.de>
//                 Susan Hert      <hert@mpi-sb.mpg.de>
//                 Lutz Kettner    <kettner@mpi-sb.mpg.de>

#ifndef CGAL_SNC_IO_PARSER_H
#define CGAL_SNC_IO_PARSER_H

#include <CGAL/license/Nef_3.h>


#include <CGAL/basic.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Nef_S2/SM_decorator.h>
#include <CGAL/Nef_3/SNC_structure.h>
#include <CGAL/Nef_3/SNC_decorator.h>
#include <CGAL/Nef_3/SNC_constructor.h>
#include <CGAL/Nef_2/Object_index.h>
#include <CGAL/Nef_S2/Normalizing.h>
#include <vector>
#include <CGAL/Fraction_traits.h>
#include <CGAL/Homogeneous_converter.h>
#include <CGAL/Cartesian_converter.h>


#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 293
#include <CGAL/Nef_2/debug.h>

#ifndef CGAL_I_DO_WANT_TO_USE_GENINFO
#include <boost/any.hpp>
#endif

#include <boost/mpl/has_xxx.hpp>

namespace CGAL {

template <class NT> class Lazy_exact_nt;
class Homogeneous_tag;
class Cartesian_tag;
template <class T> struct Simple_cartesian;
template <class T1, class T2> struct Simple_homogeneous;
template <class RT> class Quotient;

namespace Nef_3_internal{

BOOST_MPL_HAS_XXX_TRAIT_NAMED_DEF(Has_nested_Exact_kernel,Exact_kernel,false)

template <class R,bool has_exact_kernel=Has_nested_Exact_kernel<R>::value &&
                                        !std::is_floating_point<typename R::RT>::value, class FT = typename R::RT, class Kernel_tag=typename R::Kernel_tag>
struct Type_converter{
  typedef const CGAL::Point_3<R>& Point_3;
  typedef const CGAL::Vector_3<R>& Vector_3;
  typedef const CGAL::Plane_3<R>& Plane_3;

  static Point_3 convert(Point_3 p){return p;}
  static Plane_3 convert(Plane_3 p){return p;}
  static Vector_3 convert(Vector_3 v){return v;}
};

template <class R>
struct Type_converter<R, true>{
  typedef CGAL::Point_3<typename R::Exact_kernel> Point_3;
  typedef CGAL::Plane_3<typename R::Exact_kernel> Plane_3;
  typedef CGAL::Vector_3<typename R::Exact_kernel> Vector_3;

  static Point_3 convert(const CGAL::Point_3<R>& p){return p.exact();}
  static Plane_3 convert(const CGAL::Plane_3<R>& p){return p.exact();}
  static Vector_3 convert(const CGAL::Vector_3<R>& v){return v.exact();}
};

template <class R, class NT>
struct Type_converter<R, false, ::CGAL::Lazy_exact_nt<NT>, ::CGAL::Cartesian_tag >{
  typedef CGAL::Simple_cartesian< NT > EK;
  typedef CGAL::Cartesian_converter<R, EK> Converter;
  typedef CGAL::Point_3<EK> Point_3;
  typedef CGAL::Plane_3<EK> Plane_3;
  typedef CGAL::Vector_3<EK> Vector_3;

  static Point_3 convert(const CGAL::Point_3<R>& p){return Converter()(p);}
  static Plane_3 convert(const CGAL::Plane_3<R>& p){return Converter()(p);}
  static Vector_3 convert(const CGAL::Vector_3<R>& v){return Converter()(v);}
};

template <class R, class NT>
struct Type_converter<R, false, ::CGAL::Lazy_exact_nt<NT>, ::CGAL::Homogeneous_tag>{
  typedef CGAL::Simple_homogeneous< NT, CGAL::Quotient<NT> > EK;
  typedef CGAL::Homogeneous_converter<R, EK> Converter;
  typedef CGAL::Point_3<EK> Point_3;
  typedef CGAL::Plane_3<EK> Plane_3;
  typedef CGAL::Vector_3<EK> Vector_3;

  static Point_3 convert(const CGAL::Point_3<R>& p){return Converter()(p);}
  static Plane_3 convert(const CGAL::Plane_3<R>& p){return Converter()(p);}
  static Vector_3 convert(const CGAL::Vector_3<R>& v){return Converter()(v);}
};


template <class R>
typename Type_converter<R>::Point_3 get_point(const CGAL::Point_3<R>& p){ return Type_converter<R>::convert(p); }

template <class R>
typename Type_converter<R>::Plane_3 get_plane(const CGAL::Plane_3<R>& p){ return Type_converter<R>::convert(p); }

template <class R>
typename Type_converter<R>::Vector_3 get_vector(const CGAL::Vector_3<R>& v){ return Type_converter<R>::convert(v); }

} //end of Nef_3_internal

template<typename T>
class moreLeft : public T {

  typedef typename T::SHalfedge_handle  SHalfedge_handle;
  typedef typename T::Vector_3          Vector_3;
  typedef typename T::FT                FT;
  typedef typename T::RT                RT;

 public:
  moreLeft(const T& D) : T(D) {}

  int operator()(SHalfedge_handle se1, SHalfedge_handle se2) const {

    CGAL_assertion(se1 != SHalfedge_handle());
    if(se2 == SHalfedge_handle())
      return -1;

    Vector_3 vec1 = se1->circle().orthogonal_vector();
    Vector_3 vec2 = se2->circle().orthogonal_vector();

    if(vec1 == vec2)
      return 0;

    if(vec1.x() == RT(0) && vec2.x() == RT(0)) {
      if(vec1.y() != vec2.y()) {
        if(vec1.y() < vec2.y())
          return -1;
        else
          return 1;
      }
      if(vec1.z() < vec2.z())
        return -1;
      else
        return 1;
    }

    Vector_3 minus(-1,0,0);
    FT sk1(minus*vec1),  sk2(minus*vec2);
    if((sk1 >= FT(0) && sk2 <= FT(0)) ||
       (sk1 <= FT(0) && sk2 >= FT(0))) {
      if(sk1 > FT(0) || sk2 < FT(0))
        return -1;
      else
        return 1;
    }

    FT len1 = vec1.x()*vec1.x()+vec1.y()*vec1.y()+vec1.z()*vec1.z();
    FT len2 = vec2.x()*vec2.x()+vec2.y()*vec2.y()+vec2.z()*vec2.z();
    FT diff = len1*sk2*sk2 - len2*sk1*sk1;

    if(diff != FT(0)) {
      if((sk1>FT(0) && diff<FT(0)) || (sk1<FT(0) && diff>FT(0)))
        return -1;
      else
        return 1;
    }

    return 0;
  }
};

template <typename T>
class sort_vertices : public SNC_const_decorator<T> {

  typedef T SNC_structure;
  typedef CGAL::SNC_const_decorator<T>    Base;
  typedef typename T::Vertex_handle Vertex_handle;
  typedef typename T::Point_3       Point_3;

 public:
  sort_vertices(const T& D) : Base(D) {}

  bool operator() (Vertex_handle v1, Vertex_handle v2) const {
    return lexicographically_xyz_smaller(v1->point(), v2->point());
  }
};

template <typename T>
class sort_edges : public SNC_const_decorator<T> {

  typedef T SNC_structure;
  typedef CGAL::SNC_const_decorator<T>      Base;
  typedef typename T::Halfedge_handle Halfedge_handle;

 public:
  sort_edges(const T& D) : Base(D) {}

  bool operator() (Halfedge_handle e1, Halfedge_handle e2) const {
    sort_vertices<T> SORT(*this->sncp());
    if(e1->source() != e2->source())
      return SORT(e1->source(),e2->source());
    return SORT(e1->twin()->source(), e2->twin()->source());
  }
};

template <typename T>
class sort_facets : public SNC_const_decorator<T> {

  typedef T SNC_structure;
  typedef SNC_const_decorator<T>       Base;
  typedef typename T::Halffacet_handle Halffacet_handle;
  typedef typename T::SHalfedge_handle SHalfedge_handle;
  typedef typename T::Vector_3         Vector_3;
  typedef typename T::Plane_3          Plane_3;

 public:
  sort_facets(const T& D) : Base(D) {}

  bool operator() (Halffacet_handle f1, Halffacet_handle f2) const {


    Plane_3 p1(f1->plane());
    Plane_3 p2(f2->plane());

    if(p1.d() != p2.d())
      return p1.d() < p2.d();
    else if(p1.a() != p2.a())
      return p1.a() < p2.a();
    else if(p1.b() != p2.b())
      return p1.b() < p2.b();
    else if(p1.c() != p2.c())
      return p1.c() < p2.c();

    SHalfedge_handle se1 = SHalfedge_handle(f1->facet_cycles_begin());
    SHalfedge_handle se2 = SHalfedge_handle(f2->facet_cycles_begin());

    sort_vertices<T> SORT(*this->sncp());
    if(se1->source()->source() != se2->source()->source())
      return SORT(se1->source()->source(), se2->source()->source());

    se1 = se1->next();
    se2 = se2->next();
    CGAL_assertion(se1->source()->source() != se2->source()->source());
    return SORT(se1->source()->source(), se2->source()->source());


  }
};

template <typename T>
class sort_sedges : public SNC_const_decorator<T> {

  typedef T SNC_structure;
  typedef CGAL::SNC_const_decorator<T>       Base;
  typedef CGAL::SM_decorator<T>          SM_decorator;
  typedef typename T::Vertex_handle    Vertex_handle;
  typedef typename T::SHalfedge_handle SHalfedge_handle;
  typedef typename T::Sphere_circle    Sphere_circle;

 public:
  sort_sedges(const T& D) : Base(D) {}

  bool operator() (SHalfedge_handle se1, SHalfedge_handle se2) const {
    CGAL_NEF_TRACEN("sort sedges");
    if(se1 == se2) return false;
    sort_vertices<T> SORT(*this->sncp());
    CGAL_NEF_TRACEN("  center verices: " << se1->source()->source()->point() <<
                    " , " << se2->source()->source()->point());
    if(se1->source()->source() != se2->source()->source())
      return SORT(se1->source()->source(),se2->source()->source());
    if(se1 == se2->twin()) {
      if(se1->source() == se2->source()) {
        Sphere_circle vec1 = se1->circle();
        Sphere_circle vec2 = se2->circle();
        if(vec1.a() != vec2.a())
          return vec1.a() < vec2.a();
        else if(vec1.b() != vec2.b())
          return vec1.b() < vec2.b();
        return vec1.c() < vec2.c();
      }
      else
        return SORT(se1->source()->twin()->source(), se2->source()->twin()->source());
    }

    if(SORT(se1->twin()->source()->twin()->source(),
            se1->source()->twin()->source()))
      se1 = se1->twin();
    if(SORT(se2->twin()->source()->twin()->source(),
            se2->source()->twin()->source()))
      se2 = se2->twin();
    CGAL_NEF_TRACEN("  ssources " << se1->source()->twin()->source()->point()
                    << " , " << se2->source()->twin()->source()->point());
    if(se1->source() != se2->source())
      return SORT(se1->source()->twin()->source(), se2->source()->twin()->source());
    CGAL_NEF_TRACEN("  starget " << se1->twin()->source()->twin()->source()->point() <<
                    " , " << se2->twin()->source()->twin()->source()->point());
    if(se1->twin()->source()->twin()->source() != se2->twin()->source()->twin()->source())
      return SORT(se1->twin()->source()->twin()->source(), se2->twin()->source()->twin()->source());

    CGAL_assertion(se1->circle() != se2->circle());
    Sphere_circle vec1 = se1->circle();
    Sphere_circle vec2 = se2->circle();

    if(vec1.a() != vec2.a())
      return vec1.a() < vec2.a();
    else if(vec1.b() != vec2.b())
      return vec1.b() < vec2.b();
    return vec1.c() < vec2.c();
  }
};


template <typename T>
class sort_sloops : public SNC_const_decorator<T> {

  typedef T SNC_structure;
  typedef CGAL::SNC_const_decorator<T>       Base;
  typedef typename T::SHalfloop_handle SHalfloop_handle;

 public:
  sort_sloops(const T& D) : Base(D) {}

  bool operator() (SHalfloop_handle sl1, SHalfloop_handle sl2) const {
    if(sl1 == sl2) return false;
    sort_vertices<T> SORTV(*this->sncp());
    sort_facets<T> SORTF(*this->sncp());
    if(sl1->incident_sface()->center_vertex() != sl2->incident_sface()->center_vertex())
      return SORTV(sl1->incident_sface()->center_vertex(),sl2->incident_sface()->center_vertex());
    return SORTF(sl1->facet(), sl2->facet());
  }
};

template <typename T>
class sort_sface_cycle_entries : public SNC_const_decorator<T> {

  typedef T                             SNC_structure;
  typedef CGAL::SNC_const_decorator<T>  Base;
  typedef typename T::SM_decorator      SM_decorator;
  typedef typename T::Object_handle     Object_handle;
  typedef typename T::SVertex_handle    SVertex_handle;
  typedef typename T::SHalfedge_handle  SHalfedge_handle;
  typedef typename T::SHalfloop_handle  SHalfloop_handle;
  typedef typename T::SFace_handle      SFace_handle;
  typedef typename T::Point_3           Point_3;
  typedef typename T::Vector_3          Vector_3;

 public:
  sort_sface_cycle_entries(const T& D) : Base(D) {}

  bool operator() (Object_handle o1, Object_handle o2) const {
    CGAL_NEF_TRACEN("sort sface cycles ");
    SVertex_handle sv1, sv2;
    SHalfedge_handle se1, se2;
    SHalfloop_handle sl1, sl2;

    if(!CGAL::assign(se1,o1) && !CGAL::assign(sl1,o1) && !CGAL::assign(sv1,o1))
      CGAL_error_msg("wrong handle");

    if(!CGAL::assign(se2,o2) && !CGAL::assign(sl2,o2) && !CGAL::assign(sv2,o2))
      CGAL_error_msg("wrong handle");

    if(se1 != SHalfedge_handle() && se2 == SHalfedge_handle())
      return true;

    if(se1 == SHalfedge_handle() && se2 != SHalfedge_handle())
      return false;

    if(sl1 != SHalfloop_handle() && sv2 != SVertex_handle())
      return true;

    if(sl2 != SHalfloop_handle() && sv1 != SVertex_handle())
      return false;

    if(se1 != SHalfedge_handle() && se2 != SHalfedge_handle()) {
      CGAL_NEF_TRACEN("  sedges " << &*se1 << " , " << &*se2);
      sort_sedges<SNC_structure> SORT(*this->sncp());
      return SORT(se1,se2);
      /*
        sort_vertices<SNC_structure> SORT(*this->sncp());
      if(ssource(se1) != ssource(se2))
        return SORT(se1->source()->twin()->source(), se2->source()->twin()->source());
      else
        return SORT(se1->target(), se2->target());
      */
    }

    if(sl1 != SHalfloop_handle() && sl2 != SHalfloop_handle()) {
      Vector_3 vec1(sl1->circle().orthogonal_vector());
      Vector_3 vec2(sl2->circle().orthogonal_vector());
      //      CGAL_assertion(vec1 == vec2.antipode());
      if(vec1.x() != vec2.x())
        return vec1.x() < vec2.x();
      else if(vec1.y() != vec2.y())
        return vec1.y() < vec2.y();
      else if(vec1.z() != vec2.z())
        return vec1.z() < vec2.z();
    }

    CGAL_assertion(sv1 != SVertex_handle() && sv2 != SVertex_handle());
    sort_vertices<SNC_structure> SORT(*this->sncp());
    return SORT(sv1->twin()->source(), sv2->twin()->source());
  }
};

template <typename T>
class sort_sfaces : public SNC_const_decorator<T> {

  typedef T SNC_structure;
  typedef CGAL::SNC_const_decorator<T>      Base;
  typedef typename T::SM_decorator          SM_decorator;
  typedef typename T::Point_3               Point_3;
  typedef typename T::Vector_3              Vector_3;
  typedef typename T::SVertex_handle        SVertex_handle;
  typedef typename T::SHalfedge_handle      SHalfedge_handle;
  typedef typename T::SHalfloop_handle      SHalfloop_handle;
  typedef typename T::SFace_handle          SFace_handle;
  typedef typename T::SFace_cycle_iterator  SFace_cycle_iterator;
  typedef typename T::SHalfedge_around_sface_circulator
                      SHalfedge_around_sface_circulator;

 public:
  sort_sfaces(const T& D) : Base(D) {}

  bool operator() (SFace_handle sf1, SFace_handle sf2) const {
    CGAL_NEF_TRACEN("sort sfaces");
    if(&*sf1 == &*sf2) return false;

    sort_vertices<T> SORT(*this->sncp());

    CGAL_NEF_TRACEN("  vertices " << sf1->center_vertex()->point() << " , " << sf2->center_vertex()->point());
    if(sf1->center_vertex() != sf2->center_vertex())
      return SORT(sf1->center_vertex(), sf2->center_vertex());

    //    sort_sface_cycle_entries<Base> sort_cycles(*this);
    //    return sort_cycles(*sf1->sface_cycles_begin(), *sf2->sface_cycles_begin());

    SM_decorator SD(&*sf1->center_vertex());
    moreLeft<Base> ml(*this);
    Vector_3 plus(1,0,0);

    SFace_cycle_iterator fc;

    CGAL_NEF_TRACEN("  sface 1");

    SHalfedge_handle se1;
    SHalfloop_handle sl1;
    CGAL_forall_sface_cycles_of(fc,sf1) {
      if(fc.is_shalfedge()) {
        SHalfedge_handle se(fc);
        SHalfedge_around_sface_circulator ec(se),ee(se);
        CGAL_For_all(ec,ee) {
          CGAL_NEF_TRACEN("     " << ec->source()->point() <<
                 " | " << ec->circle().orthogonal_vector());
          if(ml(ec, se1) == -1)
            se1 = ec;
        }
      }
      else if(fc.is_shalfloop())
        sl1 = SHalfloop_handle(fc);
      else
        CGAL_assertion(fc.is_svertex());
    }

    CGAL_NEF_TRACEN("  sface 2");

    SHalfedge_handle se2;
    SHalfloop_handle sl2;
    CGAL_forall_sface_cycles_of(fc,sf2) {
      if(fc.is_shalfedge()) {
        SHalfedge_handle se(fc);
        SHalfedge_around_sface_circulator ec(se),ee(se);
        CGAL_For_all(ec,ee) {
          CGAL_NEF_TRACEN("     " << ec->source()->point() <<
                 " | " << ec->circle().orthogonal_vector());
          if(ml(ec, se2) == -1)
            se2 = ec;
        }
      }
      else if(fc.is_shalfloop())
        sl2 = SHalfloop_handle(fc);
      else
        CGAL_assertion(fc.is_svertex());
    }

    CGAL_NEF_TRACEN("  sedge cycles existing? " << (se1 != SHalfedge_handle())
           << " , " << (se2 != SHalfedge_handle()));

    if(se1 != SHalfedge_handle() && se2 == SHalfedge_handle())
      return true;
    if(se1 == SHalfedge_handle() && se2 != SHalfedge_handle())
      return false;

    if(se1 == SHalfedge_handle() && se2 == SHalfedge_handle()) {
      Vector_3 vec1 = sl1->circle().orthogonal_vector();
      Vector_3 vec2 = sl2->circle().orthogonal_vector();
      CGAL_NEF_TRACEN("  sloops " << vec1 << " , " << vec2);
      if(vec1.x() != vec2.x())
        return vec1.x() < vec2.x();
      else if(vec1.y() != vec2.y())
        return vec1.y() < vec2.y();
      else if(vec1.z() != vec2.z())
        return vec1.z() < vec2.z();
    }

    CGAL_assertion(se1 != SHalfedge_handle() && se2 != SHalfedge_handle());

    CGAL_NEF_TRACEN("  minimal sedge in sface 1:" << se1->source()->point() <<
           " , " << se1->circle().orthogonal_vector());
    CGAL_NEF_TRACEN("  minimal sedge in sface 2:" << se2->source()->point() <<
           " , " << se2->circle().orthogonal_vector());
    CGAL_NEF_TRACEN("result " << ml(se1,se2));
    switch(ml(se1, se2)) {
    case -1: return true;
    case  1: return false;
    }
    sort_sface_cycle_entries<T> SORTSFC(*this->sncp());
    return SORTSFC(*sf1->sface_cycles_begin(), *sf2->sface_cycles_begin());
  }
};

template <typename T>
class sort_volumes : public SNC_const_decorator<T> {

  typedef T SNC_structure;
  typedef CGAL::SNC_const_decorator<T> Base;
  typedef typename T::Volume_handle Volume_handle;
  typedef typename T::SFace_handle  SFace_handle;

 public:
  sort_volumes(const T& D) : Base(D) {}

  bool operator() (Volume_handle c1, Volume_handle c2) const {
    CGAL_NEF_TRACEN("sort volumes");
    SFace_handle sf1 = SFace_handle(c1->shells_begin());
    SFace_handle sf2 = SFace_handle(c2->shells_begin());

    sort_sfaces<T> SORT(*this->sncp());
    return SORT(sf1, sf2);
  }
};

template <typename T>
class sort_facet_cycle_entries : public T {

  typedef typename T::SNC_structure     SNC_structure;
  typedef typename T::SM_decorator      SM_decorator;
  typedef typename T::Object_handle     Object_handle;
  typedef typename T::SHalfedge_handle  SHalfedge_handle;
  typedef typename T::SHalfloop_handle  SHalfloop_handle;
  typedef typename T::SFace_handle      SFace_handle;
  typedef typename T::Point_3           Point_3;
  typedef typename T::Vector_3          Vector_3;

 public:
  sort_facet_cycle_entries(const T& D) : T(D) {}

  bool operator() (Object_handle o1, Object_handle o2) const {

    SHalfedge_handle se1, se2;
    SHalfloop_handle sl1, sl2;

    if(!CGAL::assign(se1,o1) && !CGAL::assign(sl1,o1))
      CGAL_error_msg("wrong handle");

    if(!CGAL::assign(se2,o2) && !CGAL::assign(sl2,o2))
      CGAL_error_msg("wrong handle");

    if(se1 != SHalfedge_handle() && se2 != SHalfedge_handle()) {
      sort_vertices<SNC_structure> SORT(*this->sncp());
      return SORT(se1->source()->source(), se2->source()->source());
    }

    if(se1 != SHalfedge_handle())
      return true;
    if(se2 != SHalfedge_handle())
      return false;

    CGAL_assertion(sl1 != SHalfloop_handle() &&
                        sl2 != SHalfloop_handle());

    SM_decorator SD(&*sl1->incident_sface()->center_vertex());
    Vector_3 vec1(sl1->circle().orthogonal_vector());
    Vector_3 vec2(sl2->circle().orthogonal_vector());
    //    CGAL_assertion(vec1 == vec2.antipode());
    if(vec1.x() != vec2.x())
      return vec1.x() < vec2.x();
    else if(vec1.y() != vec2.y())
      return vec1.y() < vec2.y();
    else
      return vec1.z() < vec2.z();
  }
};

template <typename T>
class sort_shell_entries : public T {

  typedef typename T::Object_handle Object_handle;
  typedef typename T::Shell_entry_iterator  Shell_entry_iterator;
  typedef typename T::SFace_handle  SFace_handle;
  typedef typename T::Point_3       Point_3;

 public:
  sort_shell_entries(const T& D) : T(D) {}

  bool operator() (Object_handle o1, Object_handle o2) const {
    SFace_handle sf1, sf2;
    CGAL::assign(sf1, o1);
    CGAL::assign(sf2, o2);
    Point_3 p1(sf1->center_vertex()->point()), p2(sf2->center_vertex()->point());
    if(p1.x() != p2.x())
      return p1.x() < p2.x();
    else if(p1.y() != p2.y())
      return p1.y() < p2.y();
    return p1.z() < p2.z();
  }
};

template<typename T>
struct find_minimal_sface_of_shell : public SNC_const_decorator<T> {

  typedef T                               SNC_structure;
  typedef CGAL::SNC_const_decorator<T>    Base;
  typedef typename T::Vertex_handle       Vertex_handle;
  typedef typename T::Halfedge_handle     Halfedge_handle;
  typedef typename T::Halffacet_handle    Halffacet_handle;
  typedef typename T::SFace_handle        SFace_handle;
  typedef typename T::SHalfedge_handle    SHalfedge_handle;
  typedef typename T::SHalfloop_handle    SHalfloop_handle;
  typedef CGAL::Unique_hash_map<SFace_handle,bool> SFace_visited_hash;

  SFace_visited_hash& Done;
  SFace_handle sf_min;
  sort_sfaces<T> SORT;

  find_minimal_sface_of_shell(const T& D, SFace_visited_hash& Vi)
    : Base(D), Done(Vi), SORT(D) {}

  void visit(SFace_handle h) {
    Done[h]=true;
    if(sf_min == SFace_handle())
      sf_min = h;
    else {
      if(SORT(h,sf_min))
        sf_min = h;
    }
  }

  void visit(Vertex_handle ) {}
  void visit(Halfedge_handle ) {}
  void visit(Halffacet_handle ) {}
  void visit(SHalfedge_handle ) {}
  void visit(SHalfloop_handle ) {}

  SFace_handle& minimal_sface() { return sf_min; }
};


template<typename Tag, typename Kernel> class Geometry_io;

template<typename Kernel>
class Geometry_io<Cartesian_tag, Kernel> {
 public:

  template <typename EK, typename K, typename Compose_> static
  typename EK::Point_3 read_point_impl(std::istream& in, Compose_) {
    typedef Fraction_traits<typename K::FT> FracTraits;
    typename FracTraits::Type hx, hy, hz, hw;
    typename FracTraits::Numerator_type num;
    typename FracTraits::Denominator_type denom(1);
    typename FracTraits::Compose composer;
    in >> num;
    hx = composer(num, denom);
    in >> num;
    hy = composer(num, denom);
    in >> num;
    hz = composer(num, denom);
    in >> num;
    hw = composer(num, denom);
    return typename EK::Point_3(hx,hy,hz,hw);
  }

  template <typename EK, typename K, typename Compose_> static
  typename EK::Plane_3 read_plane_impl(std::istream& in, Compose_) {
    typedef Fraction_traits<typename K::FT> FracTraits;
    typename FracTraits::Type a, b, c, d;
    typename FracTraits::Numerator_type num;
    typename FracTraits::Denominator_type denom(1);
    typename FracTraits::Compose composer;
    in >> num;
    a = composer(num, denom);
    in >> num;
    b = composer(num, denom);
    in >> num;
    c = composer(num, denom);
    in >> num;
    d = composer(num, denom);
    return typename EK::Plane_3(a,b,c,d);
  }

  template <typename EK, typename K> static
  typename EK::Point_3 read_point_impl(std::istream& in, Null_functor) {
    typename K::FT hx, hy, hz, hw;
    in >> hx >> hy >> hz >> hw;
    return typename EK::Point_3(hx,hy,hz,hw);
  }

  template <typename EK, typename K> static
  typename EK::Plane_3 read_plane_impl(std::istream& in, Null_functor) {
    typename K::FT a, b, c, d;
    in >> a >> b >> c >> d;
    return typename EK::Plane_3(a,b,c,d);
  }

  template <typename EK, typename K> static
  typename EK::Point_3 read_point(std::istream& in) {
    return read_point_impl<EK,K>(in, typename Fraction_traits<typename K::FT>::Compose());
  }

  template <typename EK, typename K> static
  typename EK::Plane_3 read_plane(std::istream& in) {
    return read_plane_impl<EK,K>(in, typename Fraction_traits<typename K::FT>::Compose());
  }

  template <typename R, typename Decompose_> static
  void print_point_impl(std::ostream& out, const CGAL::Point_3<R> p, Decompose_) {
    typedef Fraction_traits<typename R::FT> FracTraits;
    typedef std::vector<typename FracTraits::Numerator_type> NV;

    typename FracTraits::Numerator_type num;
    typename FracTraits::Denominator_type denom;
    typename FracTraits::Decompose decomposer;
    NV vec;

    decomposer(p.x(),num,denom);
    vec.push_back(num);
    vec.push_back(denom);
    vec.push_back(denom);
    vec.push_back(denom);
    decomposer(p.y(),num,denom);
    vec[0]*=denom;
    vec[1]*=num;
    vec[2]*=denom;
    vec[3]*=denom;
    decomposer(p.z(),num,denom);
    vec[0]*=denom;
    vec[1]*=denom;
    vec[2]*=num;
    vec[3]*=denom;
    Normalizing<Homogeneous_tag>::
      normalized(vec.begin(),vec.end());
    out << vec[0] << " " << vec[1] << " "
        << vec[2] << " " << vec[3];
  }

  template <typename R, typename Decompose_> static
  void print_vector_impl(std::ostream& out, const CGAL::Vector_3<R> p, Decompose_) {
    typedef Fraction_traits<typename R::FT> FracTraits;
    typedef typename FracTraits::Numerator_type NumType;
    typedef std::vector<NumType> NV;

    typename FracTraits::Numerator_type num;
    typename FracTraits::Denominator_type denom;
    typename FracTraits::Decompose decomposer;
    NV vec;

    decomposer(p.x(),num,denom);
    vec.push_back(num);
    vec.push_back(denom);
    vec.push_back(denom);
    decomposer(p.y(),num,denom);
    vec[0]*=denom;
    vec[1]*=num;
    vec[2]*=denom;
    decomposer(p.z(),num,denom);
    vec[0]*=denom;
    vec[1]*=denom;
    vec[2]*=num;
    Normalizing<Homogeneous_tag>::
      normalized(vec.begin(),vec.end());
    out << vec[0] << " " << vec[1] << " "
        << vec[2] << " " << NumType(1);
  }

  template <typename R, typename Decompose_> static
  void print_plane_impl(std::ostream& out, const CGAL::Plane_3<R> p, Decompose_) {
    typedef Fraction_traits<typename R::FT> FracTraits;
    typedef std::vector<typename FracTraits::Numerator_type> NV;

    typename FracTraits::Numerator_type num;
    typename FracTraits::Denominator_type denom;
    typename FracTraits::Decompose decomposer;
    NV vec;

    decomposer(p.a(),num,denom);
    vec.push_back(num);
    vec.push_back(denom);
    vec.push_back(denom);
    vec.push_back(denom);
    decomposer(p.b(),num,denom);
    vec[0]*=denom;
    vec[1]*=num;
    vec[2]*=denom;
    vec[3]*=denom;
    decomposer(p.c(),num,denom);
    vec[0]*=denom;
    vec[1]*=denom;
    vec[2]*=num;
    vec[3]*=denom;
    decomposer(p.d(),num,denom);
    vec[0]*=denom;
    vec[1]*=denom;
    vec[2]*=denom;
    vec[3]*=num;
    Normalizing<Homogeneous_tag>::
      normalized(vec.begin(),vec.end());

    out << vec[0] << " " << vec[1] << " "
        << vec[2] << " " << vec[3];
  }

  template <typename R> static
  void print_point_impl(std::ostream& out, const CGAL::Point_3<R> p, Null_functor)
  {
    out << p.x() << " " << p.y() << " " << p.z() << " " << 1;
  }

  template <typename R> static
  void print_vector_impl(std::ostream& out, const CGAL::Vector_3<R> v, Null_functor)
  {
    out << v.x() << " " << v.y() << " " << v.z() << " " << 1;
  }

  template <typename R> static
  void print_plane_impl(std::ostream& out, const CGAL::Plane_3<R> p, Null_functor)
  {
    out << p.a() << " " << p.b() << " " << p.c() << " " << p.d();
  }

  template <class R> static
  void print_point(std::ostream& out, const CGAL::Point_3<R>& p) {
    print_point_impl(out, Nef_3_internal::get_point(p), typename Fraction_traits<typename R::FT>::Decompose() );
  }

  template <class R> static
  void print_vector(std::ostream& out, const CGAL::Vector_3<R>& v) {
    print_vector_impl(out, Nef_3_internal::get_vector(v), typename Fraction_traits<typename R::FT>::Decompose() );
  }

  template <class R> static
  void print_plane(std::ostream& out, const CGAL::Plane_3<R>& p) {
    print_plane_impl(out, Nef_3_internal::get_plane(p), typename Fraction_traits<typename R::FT>::Decompose() );
  }

};

template<typename Kernel>
class Geometry_io<Homogeneous_tag, Kernel> {
 public:
  template <typename EK, typename K> static
  typename EK::Point_3 read_point(std::istream& in) {
    typename K::RT hx, hy, hz, hw;
    in >> hx >> hy >> hz >> hw;
    return typename EK::Point_3(hx, hy, hz, hw);
  }

  template <typename EK, typename K> static
  typename EK::Plane_3 read_plane(std::istream& in) {
    typename K::RT a, b, c, d;
    in >> a >> b >> c >> d;
    return typename EK::Plane_3(a, b, c, d);
  }

  template <typename R> static
  void print_point_impl(std::ostream& out, const CGAL::Point_3<R>& p) {
    out << p;
  }

  template <typename R> static
  void print_vector_impl(std::ostream& out, const CGAL::Vector_3<R>& vec) {
    out << vec;
  }

  template <typename R> static
  void print_plane_impl(std::ostream& out, const CGAL::Plane_3<R>& p) {
    out << p;
  }

  template <class R> static
  void print_point(std::ostream& out, const CGAL::Point_3<R>& p) {
    print_point_impl(out, Nef_3_internal::get_point(p) );
  }

  template <class R> static
  void print_vector(std::ostream& out, const CGAL::Vector_3<R>& v) {
    print_vector_impl(out, Nef_3_internal::get_vector(v) );
  }

  template <class R> static
  void print_plane(std::ostream& out, const CGAL::Plane_3<R>& p) {
    print_plane_impl(out, Nef_3_internal::get_plane(p) );
  }
};

template <typename SNC_structure_>
class SNC_io_parser : public SNC_decorator<SNC_structure_>
{ typedef SNC_structure_ SNC_structure;
  typedef CGAL::SNC_io_parser<SNC_structure_> Self;
  typedef CGAL::SNC_decorator<SNC_structure_> Base;
  typedef typename Base::SNC_constructor SNC_constructor;
  typedef typename SNC_structure::Sphere_map  Sphere_map;
  typedef CGAL::SM_decorator<Sphere_map>      SM_decorator;
  typedef typename SNC_structure::Infi_box    Infi_box;
  typedef typename Infi_box::Standard_kernel  Standard_kernel;
public:
  typedef typename SNC_structure::Vertex_iterator Vertex_iterator;
  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_iterator Halfedge_iterator;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::Halffacet_iterator Halffacet_iterator;
  typedef typename SNC_structure::Halffacet_handle Halffacet_handle;
  typedef typename SNC_structure::Volume_iterator Volume_iterator;
  typedef typename SNC_structure::Volume_handle Volume_handle;
  typedef typename SNC_structure::SVertex_iterator SVertex_iterator;
  typedef typename SNC_structure::SVertex_handle SVertex_handle;
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SFace_iterator SFace_iterator;
  typedef typename SNC_structure::SFace_handle SFace_handle;
  typedef typename SNC_structure::SHalfloop_iterator SHalfloop_iterator;
  typedef typename SNC_structure::SHalfloop_handle SHalfloop_handle;
  typedef typename SNC_structure::Object_iterator Object_iterator;
  typedef typename SNC_structure::Object_handle Object_handle;
  typedef typename SNC_structure::SFace_cycle_iterator SFace_cycle_iterator;
  typedef typename SNC_structure::Halffacet_cycle_iterator Halffacet_cycle_iterator;
  typedef typename SNC_structure::Shell_entry_iterator Shell_entry_iterator;
  typedef typename SNC_structure::SHalfedge_around_svertex_circulator
                                  SHalfedge_around_svertex_circulator;
  typedef typename SNC_structure::SHalfedge_around_sface_circulator
                                  SHalfedge_around_sface_circulator;
  typedef typename SNC_structure::SHalfedge_around_facet_circulator
                                  SHalfedge_around_facet_circulator;
  typedef typename SNC_structure::Point_3 Point_3;
  typedef typename SNC_structure::Plane_3 Plane_3;
  typedef typename SNC_structure::Vector_3 Vector_3;
  typedef typename SNC_structure::Sphere_point Sphere_point;
  typedef typename SNC_structure::Sphere_segment Sphere_segment;
  typedef typename SNC_structure::Sphere_circle Sphere_circle;
  typedef typename SNC_structure::Mark Mark;
  typedef typename SNC_structure::Kernel Kernel;
  typedef typename Kernel::RT RT;
  typedef typename Infi_box::Standard_point  Standard_point;
  typedef typename Infi_box::Standard_vector Standard_vector;
  typedef typename Infi_box::Standard_plane  Standard_plane;

  #ifdef CGAL_I_DO_WANT_TO_USE_GENINFO
  typedef void* GenPtr;
  #else
  typedef boost::any GenPtr;
  #endif

  using Base::visit_shell_objects;

 private:
  std::istream& in; std::ostream& out;
  bool verbose;
  bool reduce;
  bool sorted;
  bool addInfiBox;

  CGAL::Object_index<Vertex_iterator> VI;
  CGAL::Object_index<Halfedge_iterator> EI;
  CGAL::Object_index<Halffacet_iterator>    FI;
  CGAL::Object_index<Volume_iterator>   CI;
  CGAL::Object_index<SHalfedge_iterator> SEI;
  CGAL::Object_index<SHalfloop_iterator>   SLI;
  CGAL::Object_index<SFace_iterator>     SFI;
  std::list<Vertex_iterator> VL;
  std::list<Halfedge_iterator> EL;
  std::list<Halffacet_iterator> FL;
  std::list<Volume_iterator> CL;
  std::list<SHalfedge_iterator> SEL;
  std::list<SHalfloop_iterator> SLL;
  std::list<SFace_iterator> SFL;
  std::vector<Vertex_iterator>   Vertex_of;
  std::vector<Halfedge_iterator> Edge_of;
  std::vector<Halffacet_iterator>    Halffacet_of;
  std::vector<Volume_iterator>   Volume_of;
  std::vector<SVertex_iterator>  SVertex_of;
  std::vector<SHalfedge_iterator> SEdge_of;
  std::vector<SHalfloop_iterator> SLoop_of;
  std::vector<SFace_iterator>     SFace_of;
  std::size_t i,vn,en,fn,cn,sen,sln,sfn;

public:
  SNC_io_parser(std::istream& is, SNC_structure& W);
  SNC_io_parser(std::ostream& os, SNC_structure& W,
                bool sort=false, bool reduce_ = false);

  std::string index(Vertex_iterator v) const
  { return VI(v,verbose); }
  std::string index(Halfedge_iterator e) const
  { return EI(e,verbose); }
  std::string index(Halffacet_iterator f) const
  { return FI(f,verbose); }
  std::string index(Volume_iterator c) const
  { return CI(c,verbose); }
  std::string index(SHalfedge_iterator e) const
  { return SEI(e,verbose); }
  std::string index(SHalfloop_iterator l) const
  { return SLI(l,verbose); }
  std::string index(SFace_iterator f) const
  { return SFI(f,verbose); }
  std::string index(Object_iterator o) const
  { if( o == 0 )
      return this->string("undef");
    Vertex_iterator v;
    Halfedge_iterator e;
    Halffacet_iterator f;
    Volume_iterator c;
    SHalfedge_iterator se;
    SHalfloop_iterator sl;
    SFace_iterator sf;
    if( CGAL::assign( v, *o))
      return index(v);
    else if( CGAL::assign( e, *o))
      return index(e);
    else if( CGAL::assign( f, *o))
      return index(f);
    else if( CGAL::assign( c, *o))
      return index(c);
    else if( CGAL::assign( se, *o))
      return index(se);
    else if( CGAL::assign( sl, *o))
      return index(sl);
    else if( CGAL::assign( sf, *o))
      return index(sf);
    return this->string("unknown object");
  }

  bool check_sep(const char* sep) const;
  bool test_string(std::string s) const;
  void print_vertex(Vertex_handle) const;
  void print_edge(Halfedge_handle) const;
  void print_facet(Halffacet_handle) const;
  void print_volume(Volume_handle) const;
  void print_sedge(SHalfedge_handle) const;
  void print_sloop(SHalfloop_handle) const;
  void print_sface(SFace_handle) const;
  void print() const;
  void print_local_graph(Vertex_handle) const;

  template <typename NT> bool read_vertex(Vertex_handle);
  template <typename NT> bool read_edge(Halfedge_handle);
  template <typename NT> bool read_facet(Halffacet_handle);
  bool read_volume(Volume_handle);
  template <typename NT> bool read_svertex(SVertex_handle);
  template <typename NT> bool read_sedge(SHalfedge_handle);
  template <typename NT> bool read_sloop(SHalfloop_handle);
  bool read_sface(SFace_handle);
  void add_infi_box();
  void read();
  template <typename K> void read_items(int);

  static void dump(SNC_structure& W, std::ostream& os = std::cerr, bool sort = false)
  { Self O(os,W, sort); O.print(); }

  template <typename Iter, typename Index>
    void output_sorted_indexes(Iter begin, Iter end, const Index& i) const {
    int low = i[begin];
    int high = low;
    for(Iter it=begin; it != end; it++) {
      if(i[it] < low) low = i[it];
      if(i[it] > high) high = i[it];
    }
    out << low << " " << high << ", ";
  }

};

template <typename EW>
SNC_io_parser<EW>::SNC_io_parser(std::istream& is, SNC_structure& W) :
  Base(W), in(is), out(std::cout),
  reduce(false), sorted(false), addInfiBox(false),
  i(0), vn(0), en(0), fn(0), cn(0), sen(0), sln(0), sfn(0)
{
  W.clear();
  CGAL_assertion(W.is_empty());
  verbose = false;
}


template <typename EW>
SNC_io_parser<EW>::SNC_io_parser(std::ostream& os, SNC_structure& W,
                                 bool sort, bool reduce_) :
  Base(W), in(std::cin), out(os),
  addInfiBox(false),
  FI(W.halffacets_begin(),W.halffacets_end(),'F'),
  CI(W.volumes_begin(),W.volumes_end(),'C'),
  SEI(W.shalfedges_begin(),W.shalfedges_end(),'e'),
  SLI(W.shalfloops_begin(),W.shalfloops_end(),'l'),
  SFI(W.sfaces_begin(),W.sfaces_end(),'f'),
  i(0),
  vn(W.number_of_vertices()),
  en(W.number_of_halfedges()),
  fn(W.number_of_halffacets()),
  cn(W.number_of_volumes()),
  sen(W.number_of_shalfedges()),
  sln(W.number_of_shalfloops()),
  sfn(W.number_of_sfaces())
{
  verbose = (IO::get_mode(out) != CGAL::IO::ASCII &&
             IO::get_mode(out) != CGAL::IO::BINARY);
  sorted = sort;
  reduce = reduce_;
  reduce = reduce && this->is_extended_kernel() && this->is_bounded();
  sorted = sorted || reduce;

  Vertex_iterator vi;
  CGAL_forall_vertices(vi, *this->sncp()) {
    VL.push_back(vi);
    if(sorted) {
      vi->point() = normalized(vi->point());
      if(vi->has_shalfloop() &&
         sort_sloops<SNC_structure>(*this->sncp())(vi->shalfloop()->twin(),
                                                   vi->shalfloop()))
        vi->shalfloop() = vi->shalfloop()->twin();
    }
  }
  if(sorted) {
    VL.sort(sort_vertices<SNC_structure>(*this->sncp()));
  }
  if(reduce)
    for(int k=0; k<4; k++){
      VL.pop_front(); VL.pop_back();
    }
  int i = 0;
  typename std::list<Vertex_iterator>::iterator vl;
  for(vl = VL.begin(); vl != VL.end(); vl++)
    VI[*vl] = i++;

  SM_decorator SD;
  Halfedge_iterator ei;
  CGAL_forall_halfedges(ei, *this->sncp()) {
    EL.push_back(ei);
    if(sorted) {
      //      std::cerr << ei->point() << " | " << normalized(ei->point()) << " |";
      ei->point() = normalized(ei->point());
      //      std::cerr << ei->point() << std::endl;
      sort_sedges<SNC_structure> sortSE(*this->sncp());
      SHalfedge_handle new_outedge = ei->out_sedge();
      SHalfedge_around_svertex_circulator cb(new_outedge), ce(cb);
      CGAL_For_all(cb,ce) {
        if(cb != new_outedge && sortSE(cb,new_outedge))
          new_outedge = cb;
      }
      ei->out_sedge() = new_outedge;
    }
  }
  if(sorted) EL.sort(sort_edges<SNC_structure>(*this->sncp()));
  if(reduce)
    for(int k=0; k<12; k++){
      EL.pop_front(); EL.pop_back();
    }
  i = 0;
  typename std::list<Halfedge_iterator>::iterator el;
  for(el = EL.begin(); el != EL.end(); el++)
    EI[*el] = i++;

  Halffacet_iterator fi;
  CGAL_forall_halffacets(fi, *this->sncp()){
    if(sorted) {
      sort_sedges<SNC_structure> sortSE(*this->sncp());
      Halffacet_cycle_iterator fc;
      for(fc = fi->facet_cycles_begin();
          fc != fi->facet_cycles_end(); ++fc) {
        if(fc.is_shalfedge()) {
          SHalfedge_handle se(fc);
          if(this->sncp()->is_boundary_object(se))
            this->sncp()->undef_boundary_item(se);
          SHalfedge_around_facet_circulator sfc(fc), send(sfc);
          CGAL_For_all(sfc, send) {
            if(sortSE(sfc, se))
              se = sfc;
          }
          this->sncp()->store_boundary_item(se,fc);
          *fc = make_object(se);
        }
      }
      fi->plane() = normalized(fi->plane());
      fi->boundary_entry_objects().sort(sort_facet_cycle_entries<Base>(*this));
    }
    FL.push_back(fi);
  }
  if(sorted) FL.sort(sort_facets<SNC_structure>(*this->sncp()));
  if(reduce) {
    for(int k=0; k<6; k++){
      FL.pop_front();
      FL.pop_back();
    }
  }
  i = 0;
  typename std::list<Halffacet_iterator>::iterator fl;
  for(fl = FL.begin(); fl != FL.end(); fl++)
    FI[*fl] = i++;

  SHalfedge_iterator sei;
  CGAL_forall_shalfedges(sei, *this->sncp()) {
    SEL.push_back(sei);
    if(sorted)
      sei->circle() = normalized(sei->circle());
  }
  if(sorted) SEL.sort(sort_sedges<SNC_structure>(*this->sncp()));
  if(reduce)
    for(int k=0; k<24; k++){
      SEL.pop_front(); SEL.pop_back();
    }
  i = 0;
  typename std::list<SHalfedge_iterator>::iterator sel;
  for(sel = SEL.begin(); sel != SEL.end(); sel++)
    SEI[*sel] = i++;

  SHalfloop_iterator sli;
  CGAL_forall_shalfloops(sli, *this->sncp()) {
    SLL.push_back(sli);
    if(sorted)
      sli->circle() = normalized(sli->circle());
  }
  if(sorted) SLL.sort(sort_sloops<SNC_structure>(*this->sncp()));
  i = 0;
  typename std::list<SHalfloop_iterator>::iterator sll;
  for(sll = SLL.begin(); sll != SLL.end(); sll++)
    SLI[*sll] = i++;

  SFace_iterator sfi;
  CGAL_forall_sfaces(sfi, *this->sncp()) {
    if(sorted) {
      SFace_cycle_iterator fc;
      CGAL_forall_sface_cycles_of(fc, sfi) {
        if(fc.is_shalfedge()) {
          SHalfedge_handle se(fc);
          if(this->sncp()->is_sm_boundary_object(se))
            this->sncp()->undef_sm_boundary_item(se);
          SHalfedge_around_sface_circulator cb(se), ce(cb);
          CGAL_For_all(cb,ce) {
            if(cb->source() != se->source()) {
              if(lexicographically_xyz_smaller(cb->source()->twin()->source()->point(),
                                               se->source()->twin()->source()->point()))
                se = cb;
            }
            else
              if(lexicographically_xyz_smaller(cb->twin()->source()->twin()->source()->point(),
                                               se->twin()->source()->twin()->source()->point()))
                se = cb;
          }
          this->sncp()->store_sm_boundary_item(se,fc);
          *fc = make_object(se);
        }
      }
      sfi->boundary_entry_objects().sort(sort_sface_cycle_entries<Base>(*this));
    }
    SFL.push_back(sfi);
  }
  if(sorted) SFL.sort(sort_sfaces<SNC_structure>(*this->sncp()));
  if(reduce)
    for(int k=0; k<8; k++){
      SFL.pop_front(); SFL.pop_back();
    }
  i = 0;
  typename std::list<SFace_iterator>::iterator sfl;
  for(sfl = SFL.begin(); sfl != SFL.end(); sfl++)
    SFI[*sfl] = i++;

  Volume_iterator ci;
  CGAL::Unique_hash_map<SFace_handle,bool> Done(false);
  find_minimal_sface_of_shell<SNC_structure> findMinSF(*this->sncp(),Done);
  CGAL_forall_volumes(ci, *this->sncp()) {
    if(sorted) {
      Shell_entry_iterator it;
      CGAL_forall_shells_of(it,ci) {
        findMinSF.minimal_sface() = SFace_handle(it);
        visit_shell_objects(SFace_handle(it),findMinSF);
        *it = make_object(findMinSF.minimal_sface());
      }
      ci->shell_entry_objects().sort(sort_shell_entries<Base>(*this));
    }
    CL.push_back(ci);
  }

  if(sorted) CL.sort(sort_volumes<SNC_structure>(*this->sncp()));
  if(reduce)
    CL.pop_front();
  i = 0;
  typename std::list<Volume_iterator>::iterator cl;
  for(cl = CL.begin(); cl != CL.end(); cl++)
    CI[*cl] = i++;

  VI[W.vertices_end()]=-2;
  EI[W.halfedges_end()]=-2;
  FI[W.halffacets_end()]=-2;
  CI[W.volumes_end()]=-2;
  SEI[W.shalfedges_end()]=-2;
  SLI[W.shalfloops_end()]=-2;
  SFI[W.sfaces_end()]=-2;
}

template <typename EW>
bool SNC_io_parser<EW>::check_sep(const char* sep) const
{
  char c;
  do in.get(c); while (isspace(c));
  while (*sep != '\0') {
    if (*sep != c) {
      in.putback(c);
      return false;
    }
    ++sep; in.get(c);
  }
  in.putback(c);
  return true;
}

template <typename EW>
bool SNC_io_parser<EW>::test_string(std::string s) const {
  std::string s2;
  in >> s2;
  return (s==s2);
}

template <typename EW>
void SNC_io_parser<EW>::print() const
{
  out << "Selective Nef Complex" << std::endl;
  if(this->is_extended_kernel() && (!reduce || !this->is_bounded()))
    out << "extended" << std::endl;
  else
    out << "standard" << std::endl;
  out << "vertices   " << VL.size() << std::endl;
  out << "halfedges  " << EL.size() << std::endl;
  out << "facets     " << FL.size() << std::endl;
  out << "volumes    " << CL.size() << std::endl;
  out << "shalfedges " << SEL.size() << std::endl;
  out << "shalfloops " << SLL.size() << std::endl;
  out << "sfaces     " << SFL.size() << std::endl;

  if (verbose)
    out << "/* Vertex: index { svs sve ses see sfs sfe sl,"
        << " mark, point } */\n";
  typename std::list<Vertex_iterator>::const_iterator v;
  for(v=VL.begin();v!=VL.end();v++)
    print_vertex(*v);

  if (verbose)
  out << "/* Edge: index { twin, source, isolated incident_object,"
      << " mark } */\n";
  typename std::list<Halfedge_iterator>::const_iterator e;
  for(e=EL.begin();e!=EL.end();e++)
    print_edge(*e);

  if (verbose)
  out << "/* Facet: index { twin, fclist, ivlist, volume | plane } mark */\n";
  typename std::list<Halffacet_iterator>::const_iterator f;
  for(f=FL.begin();f!=FL.end();f++)
    print_facet(*f);

  if (verbose)
  out << "/* Volume: index { shlist } mark  */\n";
  typename std::list<Volume_iterator>::const_iterator c;
  for(c=CL.begin();c!=CL.end();c++)
    print_volume(*c);

  if (verbose)
  out << "/* SEdge: index { twin, sprev, snext, source, sface,"
      << " prev, next, facet } */\n";
  typename std::list<SHalfedge_iterator>::const_iterator se;
  for(se=SEL.begin();se!=SEL.end();se++)
    print_sedge(*se);

  if (verbose)
  out << "/* SLoop: index { twin, sface, facet } */" << std::endl;
  typename std::list<SHalfloop_iterator>::const_iterator sl;
  for(sl=SLL.begin();sl!=SLL.end();sl++)
    print_sloop(*sl);

  if (verbose)
  out << "/* SFace: index { fclist, ivlist, sloop, volume } */" << std::endl;
  typename std::list<SFace_iterator>::const_iterator sf;
  for(sf=SFL.begin();sf!=SFL.end();sf++)
    print_sface(*sf);

  out << "/* end Selective Nef complex */" << std::endl;
}

template <typename EW>
void SNC_io_parser<EW>::read()
{
  if ( !check_sep("Selective Nef Complex") )
  {
    CGAL_warning_msg(false, "SNC_io_parser::read: no SNC header.");
    return;
  }
  std::string kernel_type;
  in >> kernel_type;
  CGAL_assertion(kernel_type == "standard" || kernel_type == "extended");
  if ( !(check_sep("vertices") && (in >> vn)) )
  {
    CGAL_warning_msg(false, "SNC_io_parser::read: wrong vertex line.");
    return;
  }
  if ( !(check_sep("halfedges") && (in >> en) && (en%2==0)) )
  {
    CGAL_warning_msg(false, "SNC_io_parser::read: wrong edge line.");
    return;
  }
  if ( !(check_sep("facets") && (in >> fn) && (fn%2==0)) )
  {
    CGAL_warning_msg(false, "SNC_io_parser::read: wrong facet line.");
  }
  if ( !(check_sep("volumes") && (in >> cn)) )
  {
    CGAL_warning_msg(false, "SNC_io_parser::read: wrong volume line.");
    return;
  }
  if ( !(check_sep("shalfedges") && (in >> sen)) )
  {
   CGAL_warning_msg(false, "SNC_io_parser::read: wrong sedge line.");
   return;
  }
  if ( !(check_sep("shalfloops") && (in >> sln)) )
  {
    CGAL_warning_msg(false, "SNC_io_parser::read: wrong sloop line.");
    return;
  }
  if ( !(check_sep("sfaces") && (in >> sfn)) )
  {
    CGAL_warning_msg(false, "SNC_io_parser::read: wrong sface line.");
    return;
  }

  addInfiBox = (kernel_type == "standard" && Infi_box::extended_kernel());

  for(i=0; i<vn; ++i)  Vertex_of.push_back(this->sncp()->new_vertex_only());
  for(i=0; i<en; ++i)  Edge_of.push_back(this->sncp()->new_halfedge_only());
  for(i=0; i<fn; ++i)  Halffacet_of.push_back(this->sncp()->new_halffacet_only());
  for(i=0; i<cn; ++i)  Volume_of.push_back(this->sncp()->new_volume_only());
  for(i=0; i<sen; ++i) SEdge_of.push_back(this->sncp()->new_shalfedge_only());
  for(i=0; i<sln; ++i) SLoop_of.push_back(this->sncp()->new_shalfloop_only());
  for(i=0; i<sfn; ++i) SFace_of.push_back(this->sncp()->new_sface_only());

  if(addInfiBox) {
    Volume_of.push_back(this->sncp()->new_volume_only());
    read_items<Standard_kernel>(1);
    add_infi_box();
  } else
    read_items<Kernel>(0);
}

template <typename EW>
template <typename K>
void SNC_io_parser<EW>::read_items(int plus01) {

  typename std::vector<Vertex_iterator>::iterator vi;
  for(vi=Vertex_of.begin(); vi!=Vertex_of.end(); ++vi) {
    if (!read_vertex<K>(*vi))
    {
      CGAL_warning_msg(false, "SNC_io_parser::read: error in node line");
      return;
    }
  }

  typename std::vector<Halfedge_iterator>::iterator ei;
  for(ei=Edge_of.begin(); ei!=Edge_of.end(); ++ei) {
    if (!read_edge<K>(*ei))
    {
      CGAL_warning_msg(false, "SNC_io_parser::read: error in edge line");
      return;
    }
  }

  typedef typename std::vector<Halffacet_iterator>::iterator vhf_iterator;
  vhf_iterator fi;
  for(fi=Halffacet_of.begin(); fi!=Halffacet_of.end(); ++fi) {
    if (!read_facet<K>(*fi))
    {
      CGAL_warning_msg(false, "SNC_io_parser::read: error in facet line");
      return;
    }
  }
  typename std::vector<Volume_iterator>::iterator ci;
  for(ci=Volume_of.begin()+plus01; ci!=Volume_of.end(); ++ci) {
    if (!read_volume(*ci))
    {
      CGAL_warning_msg(false, "SNC_io_parser::read: error in volume line");
      return;
    }
  }
  typename std::vector<SHalfedge_iterator>::iterator sei;
  for(sei=SEdge_of.begin(); sei!=SEdge_of.end(); ++sei) {
    if (!read_sedge<K>(*sei))
    {
      CGAL_warning_msg(false, "SNC_io_parser::read: error in sedge line");
      return;
    }
  }
  typename std::vector<SHalfloop_iterator>::iterator sli;
  for(sli=SLoop_of.begin(); sli!=SLoop_of.end(); ++sli) {
    if (!read_sloop<K>(*sli))
    {
      CGAL_warning_msg(false, "SNC_io_parser::read: error in sloop line");
      return;
    }
  }
  typename std::vector<SFace_iterator>::iterator sfi;
  for(sfi=SFace_of.begin(); sfi!=SFace_of.end(); ++sfi) {
    if (!read_sface(*sfi))
    {
      CGAL_warning_msg(false, "SNC_io_parser::read: error in sface line");
      return;
    }
  }

  SNC_constructor C(*this->sncp());
  C.assign_indices();
}


template <typename EW>
void SNC_io_parser<EW>::print_vertex(Vertex_handle v) const
{ // syntax: index { svs sve, ses see, sfs sfe, sl | point } mark
  SM_decorator SD(&*v);
  out << index(v) << " { ";
  if(sorted) {

    output_sorted_indexes(v->svertices_begin(),
                          v->svertices_end(), EI);
    output_sorted_indexes(v->shalfedges_begin(),
                          v->shalfedges_end(), SEI);
    output_sorted_indexes(v->sfaces_begin(),
                          v->sfaces_end(), SFI);
    out << index(SD.shalfloop()) << " | ";
  }
  else {
    out
    << index(v->svertices_begin()) << " "
    << index(v->svertices_last()) << ", "
    << index(v->shalfedges_begin()) << " "
    << index(v->shalfedges_last()) << ", "
    << index(v->sfaces_begin()) << " "
    << index(v->sfaces_last()) << ", "
    << index(SD.shalfloop()) << " | ";
  }
  if(reduce) {
    Geometry_io<typename Standard_kernel::Kernel_tag, Standard_kernel>::
      print_point(out, Infi_box::standard_point(v->point()));
  }
  else
    Geometry_io<typename Kernel::Kernel_tag, Kernel>::print_point(out, v->point());
  out << " } "  << v->mark() << std::endl;
}

template <typename EW>
template <typename K>
bool SNC_io_parser<EW>::
read_vertex(Vertex_handle vh) {

  bool OK = true;
  int index;
  #ifdef CGAL_NEF_NATURAL_COORDINATE_INPUT
  typename K::RT hx, hy, hz, hw;
  #endif

  in >> index;
  OK = OK && test_string("{");
  vh->sncp() = this->sncp();

  in >> index;
  if(index >= (int)en)
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  vh->svertices_begin() = (index >= 0 ? Edge_of[index] : this->svertices_end());
  in >> index;
  if(index >= int(en))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  vh->svertices_last()  = index >= 0 ? Edge_of[index] : this->svertices_end();
  OK = OK && test_string(",");
  in >> index;
  if(index >= int(sen))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  vh->shalfedges_begin() = index >= 0 ? SEdge_of[index] : this->shalfedges_end();
  in >> index;
  if(index >= int(sen))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  vh->shalfedges_last()  = index >= 0 ? SEdge_of[index] : this->shalfedges_end();
  OK = OK && test_string(",");
  in >> index;
  if(index >= int(sfn))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  vh->sfaces_begin() = index >= 0 ? SFace_of[index] : this->sfaces_end();
  in >> index;
  if(index >= int(sfn))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  vh->sfaces_last()  = index >= 0 ? SFace_of[index] : this->sfaces_end();
  OK = OK && test_string(",");
  in >> index;
  if(index >= int(sln))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  vh->shalfloop() = index >= 0 ? SLoop_of[index] : this->shalfloops_end();
  OK = OK && test_string("|");
#ifdef CGAL_NEF_NATURAL_COORDINATE_INPUT
  in >> hx >> hy >> hz >> hw;
  vh->point() = Point_3(hx,hy,hz,hw);
#else
  vh->point() =
    Geometry_io<typename K::Kernel_tag, Kernel>::template read_point<Kernel, K>(in);
#endif
  OK = OK && test_string("}");
  in >> vh->mark();

  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::print_edge(Halfedge_handle e) const
{ // syntax: index { twin, source, isolated incident_object | spoint } mark
  SM_decorator D(&*e->source());
  out << index(e) << " { " << index(e->twin()) << ", "
      << index(e->source()) << ", ";
  if ( D.is_isolated(e) ) out << "1 " << index(e->incident_sface());
  else out << "0 " << index(e->out_sedge());
  out << " | ";
  if(reduce) {
    Standard_point sp = Infi_box::standard_point(e->point());
    Geometry_io<typename Standard_kernel::Kernel_tag, Standard_kernel>::
      print_vector(out, sp-CGAL::ORIGIN);
  }
  else
    Geometry_io<typename Kernel::Kernel_tag, Kernel>::
      print_vector(out, e->vector());

  out << " } "<< e->mark();
#ifdef CGAL_NEF_OUTPUT_INDEXES
  out << " " << e->get_index();
#endif
  out << std::endl;
}

template <typename EW>
template <typename K>
bool SNC_io_parser<EW>::
read_edge(Halfedge_handle eh) {

  bool OK = true;
  int index;
#ifdef CGAL_NEF_NATURAL_COORDINATE_INPUT
  typename K::RT hx,hy,hz,hw;
#endif
  in >> index;
  OK = OK && test_string("{");

  in >> index;
  if(index < 0 || index >= int(en))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  eh->twin() = Edge_of[index];
  OK = OK && test_string(",");
  in >> index;
  if(index < 0 || index >= int(vn))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  eh->center_vertex() = Vertex_of[index];
  OK = OK && test_string(",");
  in >> index;
  if(index == 0) {
    in >> index;
    if(index < 0 || index >= int(sen))
    {
      in.setstate(std::ios_base::badbit);
      return false;
    }
    eh->out_sedge() = SEdge_of[index];
  } else {
    in >> index;
    if(index < 0 || index >= int(sfn))
    {
      in.setstate(std::ios_base::badbit);
      return false;
    }
    eh->incident_sface() = SFace_of[index];
  }
  OK = OK && test_string("|");
#ifdef CGAL_NEF_NATURAL_COORDINATE_INPUT
  in >> hx >> hy >> hz >> hw;
  eh->point() = Sphere_point(hx,hy,hz);
#else
  eh->point() =
    Geometry_io<typename K::Kernel_tag, Kernel>::template read_point<Kernel,K>(in);
#endif
  OK = OK && test_string("}");
  in >> eh->mark();

  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::print_facet(Halffacet_handle f) const
{ // syntax: index { twin, fclist, ivlist, volume | plane } mark
  out << index(f) << " { ";
  out << index(f->twin()) << ", ";
  Halffacet_cycle_iterator it;
  CGAL_forall_facet_cycles_of(it,f)
    if ( it.is_shalfedge() ) out << index(SHalfedge_handle(it)) << ' ';
  out << ", ";
  CGAL_forall_facet_cycles_of(it,f)
    if ( it.is_shalfloop() ) out << index(SHalfloop_handle(it)) << ' ';
  out << ", " << index(f->incident_volume()) << " | ";
  if(reduce) {
    Geometry_io<typename Standard_kernel::Kernel_tag, Standard_kernel>::
      print_plane(out, Infi_box::standard_plane(f->plane()));
  }
  else
    Geometry_io<typename Kernel::Kernel_tag, Kernel>::print_plane(out, f->plane());

  out << " } " << f->mark() << std::endl;
}

template <typename EW>
template <typename K>
bool SNC_io_parser<EW>::
read_facet(Halffacet_handle fh) {

  bool OK = true;
  int index;
  char cc;
#ifdef CGAL_NEF_NATURAL_COORDINATE_INPUT
  typename K::RT a,b,c,d;
#endif

  in >> index;
  OK = OK && test_string("{");

  in >> index;
  if(index < 0 || index >= int(fn))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  fh->twin() = Halffacet_of[index];
  OK = OK && test_string(",");

  in >> cc;
  while(isdigit(cc)) {
    in.putback(cc);
    in >> index;
    if(index < 0 || index >= int(sen))
    {
      in.setstate(std::ios_base::badbit);
      return false;
    }
    fh->boundary_entry_objects().push_back(make_object(SEdge_of[index]));
    in >> cc;
  }

  in >> cc;
  while(isdigit(cc)) {
    in.putback(cc);
    in >> index;
    if(index < 0 || index >= int(sln))
    {
      in.setstate(std::ios_base::badbit);
      return false;
    }
    fh->boundary_entry_objects().push_back(make_object(SLoop_of[index]));
    in >> cc;
  }

  in >> index;
  if(index < 0  || index >= int(cn))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  fh->incident_volume() = Volume_of[index+addInfiBox];
  OK = OK && test_string("|");
#ifdef CGAL_NEF_NATURAL_COORDINATE_INPUT
  in >> a >> b >> c >> d;
  fh->plane() = Plane_3(a,b,c,d);
#else
  fh->plane() =
    Geometry_io<typename K::Kernel_tag, Kernel>::
    template read_plane<Kernel, K>(in);
#endif
  OK = OK && test_string("}");
  in >> fh->mark();

  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::print_volume(Volume_handle c) const
{ // syntax: index { shlist } mark
  out << index(c) << " { ";
  Shell_entry_iterator it;
  CGAL_forall_shells_of(it,c)
    if(!reduce || Infi_box::is_standard(SFace_handle(it)->center_vertex()->point()))
      out << index(SFace_handle(it)) << ' ';
  out << "} " << c->mark() << std::endl;
}

template <typename EW>
bool SNC_io_parser<EW>::
read_volume(Volume_handle ch) {

  bool OK = true;
  int index;
  char cc;

  in >> index;
  OK = OK && test_string("{");

  in >> cc;
  while(isdigit(cc)) {
    in.putback(cc);
    in >> index;
    if(index < 0 || index >= int(sfn))
    {
      in.setstate(std::ios_base::badbit);
      return false;
    }
    ch->shell_entry_objects().push_back(make_object(SFace_of[index]));
    in >> cc;
  }
  in >> ch->mark();

  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::
print_sedge(SHalfedge_handle e) const {
//index { twin, sprev, snext, source, sface, prev, next, facet | circle } mark
  out << index(e) << " { "
      << index(e->twin()) << ", "
      << index(e->sprev()) << ", " << index(e->snext()) << ", "
      << index(e->source()) << ", " << index(e->incident_sface()) << ", "
      << index(e->prev()) << ", " << index(e->next()) << ", "
      << index(e->facet())
      << " | ";
  if(reduce) {
    Geometry_io<typename Standard_kernel::Kernel_tag, Standard_kernel>::
      print_plane(out, Infi_box::standard_plane(e->circle()));
  }
  else
    Geometry_io<typename Kernel::Kernel_tag, Kernel>::
      print_plane(out, (Plane_3) e->circle());

  out << " } " << e->mark();
#ifdef CGAL_NEF_OUTPUT_INDEXES
  out << " " << e->get_forward_index()
      << " " << e->get_backward_index();
#endif
  out << std::endl;
}

template <typename EW>
template <typename K>
bool SNC_io_parser<EW>::
read_sedge(SHalfedge_handle seh) {

  bool OK = true;
  int index;
#ifdef CGAL_NEF_NATURAL_COORDINATE_INPUT
  typename K::RT a,b,c,d;
#endif

  in >> index;
  OK = OK && test_string("{");

  in >> index;
  if(index < 0 || index >= int(sen))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  seh->twin() = SEdge_of[index];
  OK = OK && test_string(",");
  in >> index;
  if(index < 0 || index >= int(sen))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  seh->sprev() = SEdge_of[index];
  OK = OK && test_string(",");
  in >> index;
  if(index < 0 || index >= int(sen))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  seh->snext() = SEdge_of[index];
  OK = OK && test_string(",");
  in >> index;
  if(index < 0 || index >= int(en))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  seh->source() = Edge_of[index];
  OK = OK && test_string(",");
  in >> index;
  if(index < 0 || index >= int(sfn))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  seh->incident_sface() = SFace_of[index];
  OK = OK && test_string(",");
  in >> index;
  if(index < 0 || index >= int(sen))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  seh->prev() = SEdge_of[index];
  OK = OK && test_string(",");
  in >> index;
  if(index < 0 || index >= int(sen))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  seh->next() = SEdge_of[index];
  OK = OK && test_string(",");
  in >> index;
  if(index < 0 || index >= int(fn))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  seh->facet() = Halffacet_of[index];
  OK = OK && test_string("|");
#ifdef CGAL_NEF_NATURAL_COORDINATE_INPUT
  in >> a >> b >> c >> d;
  seh->circle() = Sphere_circle(Plane_3(a,b,c,d));
#else
  seh->circle() =
    Geometry_io<typename K::Kernel_tag, Kernel>::
    template read_plane<Kernel, K>(in);
#endif
  OK = OK && test_string("}");
  in >> seh->mark();

  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::
print_sloop(SHalfloop_handle l) const
{ // syntax: index { twin, sface, facet | circle } mark
  out << index(l) << " { "
      << index(l->twin()) << ", " << index(l->incident_sface()) << ", "
      << index(l->facet())
      << " | ";
  if(reduce) {
    Geometry_io<typename Standard_kernel::Kernel_tag, Standard_kernel>::
      print_plane(out, Infi_box::standard_plane(l->circle()));
  }
  else
    Geometry_io<typename Kernel::Kernel_tag, Kernel>::
      print_plane(out, (Plane_3) l->circle());

  out << " } " << l->mark() << "\n";
}

template <typename EW>
template <typename K>
bool SNC_io_parser<EW>::
read_sloop(SHalfloop_handle slh) {

  bool OK = true;
  int index;
#ifdef CGAL_NEF_NATURAL_COORDINATE_INPUT
  typename K::RT a,b,c,d;
#endif

  in >> index;
  OK = OK && test_string("{");

  in >> index;
  if(index < 0 || index >= (int)(sln))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  slh->twin() = SLoop_of[index];
  OK = OK && test_string(",");
  in >> index;
  if(index < 0 || index >= (int)(sfn))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  slh->incident_sface() = SFace_of[index];
  OK = OK && test_string(",");
  in >> index;
  if(index < 0 || index >= (int)(fn))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  slh->facet() = Halffacet_of[index];
  OK = OK && test_string("|");
#ifdef CGAL_NEF_NATURAL_COORDINATE_INPUT
  in >> a >> b >> c >> d;
  slh->circle() = Sphere_circle(Plane_3(a,b,c,d));
#else
  slh->circle() =
    Geometry_io<typename K::Kernel_tag, Kernel>::
    template read_plane<Kernel, K>(in);
#endif
  OK = OK && test_string("}");
  in >> slh->mark();

  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::
print_sface(SFace_handle f) const
{ // syntax: index { vertex, fclist, ivlist, sloop, volume }
  SM_decorator D(&*f->center_vertex());
  out << index(f) << " { " << index(f->center_vertex()) << ", ";
  SFace_cycle_iterator it;
  CGAL_forall_sface_cycles_of(it,f)
    if ( it.is_shalfedge() ) out << index(SHalfedge_handle(it)) << ' ';
  out << ", ";
  CGAL_forall_sface_cycles_of(it,f)
    if ( it.is_svertex() ) out << index(SVertex_handle(it)) << ' ';
  out << ", ";
  CGAL_forall_sface_cycles_of(it,f)
    if ( it.is_shalfloop() ) out << index(SHalfloop_handle(it));
  out << ", " << index(f->volume()) << " } " << f->mark() <<"\n";
}

template <typename EW>
bool SNC_io_parser<EW>::
read_sface(SFace_handle sfh) {

  bool OK = true;
  int index;
  char cc;

  in >> index;
  OK = OK && test_string("{");

  in >> index;
  if(index < 0 || index >= (int)(vn))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  sfh->center_vertex() = Vertex_of[index];
  OK = OK && test_string(",");

  in >> cc;
  while(isdigit(cc)) {
    in.putback(cc);
    in >> index;
    //    sfh->boundary_entry_objects().push_back(SEdge_of[index]);
    SM_decorator SD(&*sfh->center_vertex());
    if(index < 0 || index >= (int)(sen))
    {
      in.setstate(std::ios_base::badbit);
      return false;
    }
    SD.link_as_face_cycle(SEdge_of[index],sfh);
    in >> cc;
  }

  in >> cc;
  while(isdigit(cc)) {
    in.putback(cc);
    in >> index;
    if(index < 0 || index >= (int)(en))
    {
      in.setstate(std::ios_base::badbit);
      return false;
    }
    sfh->boundary_entry_objects().push_back(make_object(Edge_of[index]));
    this->sncp()->store_sm_boundary_item(Edge_of[index], --(sfh->sface_cycles_end()));
    in >> cc;
  }

  in >> cc;
  while(isdigit(cc)) {
    in.putback(cc);
    in >> index;
    if(index < 0 || index >= (int)(sln))
    {
      in.setstate(std::ios_base::badbit);
      return false;
    }
    sfh->boundary_entry_objects().push_back(make_object(SLoop_of[index]));
    this->sncp()->store_sm_boundary_item(SLoop_of[index], --(sfh->sface_cycles_end()));
    in >> cc;
  }

  in >> index;
  if(index < 0 || index >= (int)(cn))
  {
    in.setstate(std::ios_base::badbit);
    return false;
  }
  sfh->volume() = Volume_of[index+addInfiBox];
  OK = OK && test_string("}");
  in >> sfh->mark();

  return OK;
}

template <typename EW>
void SNC_io_parser<EW>::print_local_graph(Vertex_handle v) const
{ SM_decorator D(&*v);
  out << "Local Graph "
      << D.number_of_vertices() << " " << D.number_of_edges() << " "
      << D.number_of_loops() << " " << D.number_of_faces() << " "
      << std::endl;
  if (verbose)
    out << "/* index { twin, source, isolated incident_object, mark } */\n";
  SVertex_iterator vit;
  CGAL_forall_svertices_of(vit,v) print_edge(vit);
  if (verbose)
    out << "/* index { twin, sprev, snext, source, sface,"
        << " prev, next, facet } */\n";
  SHalfedge_iterator eit;
  CGAL_forall_shalfedges_of(eit,v) print_sedge(eit);
  if (verbose)
    out << "/* index { twin, sface, facet } */" << std::endl;
  if ( D.has_sloop() )
  { print_sloop(D.loop()); print_sloop(twin(D.loop())); }
  if (verbose)
    out << "/* index { fclist, ivlist, sloop, volume } */" << std::endl;
  SFace_iterator fit;
  CGAL_forall_sfaces_of(fit,v) print_sface(fit);
  out.flush();
}

template <typename EW>
void SNC_io_parser<EW>::add_infi_box() {

  for(i=0; i<8; ++i)  Vertex_of.push_back(this->sncp()->new_vertex_only());
  for(i=0; i<24; ++i)  Edge_of.push_back(this->sncp()->new_halfedge_only());
  for(i=0; i<12; ++i)  Halffacet_of.push_back(this->sncp()->new_halffacet_only());
  for(i=0; i<48; ++i) SEdge_of.push_back(this->sncp()->new_shalfedge_only());
  for(i=0; i<16; ++i) SFace_of.push_back(this->sncp()->new_sface_only());

  typename Standard_kernel::RT hx,hy,hz,hw;
  for(int i=0; i<8; ++i) {
    Vertex_handle vh = Vertex_of[vn+i];
    vh->svertices_begin() = Edge_of[en+3*i];
    vh->svertices_last()  = Edge_of[en+3*i+2];
    vh->shalfedges_begin() = SEdge_of[sen+6*i];
    vh->shalfedges_last()  = SEdge_of[sen+6*i+5];
    vh->sfaces_begin() = SFace_of[sfn+2*i];
    vh->sfaces_last()  = SFace_of[sfn+2*i+1];
    vh->shalfloop() = this->shalfloops_end();
    hx = i % 2 ? -1 : 1;
    hy = i % 4 > 1 ? -1 : 1;
    hz = i > 3 ? -1 : 1;
    vh->point() = Infi_box::create_extended_point(hx, hy, hz);
    vh->mark() = 1;
    vh->sncp() = this->sncp();
  }

  int seOff[3] = {0, 1, 3};
  int twinIdx[24] = { 3, 7,14,
                      0,10,17,
                      9, 1,20,
                      6, 4,23,
                      15,19, 2,
                      12,22, 5,
                      21,13, 8,
                      18,16,11};

  for(int i = 0; i < 24; ++i) {
    Halfedge_handle eh = Edge_of[en+i];
    eh->twin() = Edge_of[en+twinIdx[i]];
    eh->center_vertex() = Vertex_of[vn+(i/3)];
    eh->out_sedge() = SEdge_of[sen+(i/3*6)+seOff[i%3]];
    switch(i%3) {
    case 0 :
      hx = i % 6 ? 1 : -1;
      hy = hz = 0;
      break;
    case 1:
      hy = i % 12 >= 6 ? 1 : -1;
      hx = hz = 0;
      break;
    case 2:
      hz = i >= 12 ? 1 : -1;
      hx = hy = 0;
      break;
    }
    eh->point() = Sphere_point(hx,hy,hz);
    eh->mark() = 1;
  }

  int bnd[12] = {19, 18, 43, 42, 35, 34,
                 47, 46, 39, 38, 45, 44};
  for(int i = 0; i < 12; ++i) {
    Halffacet_handle fh = Halffacet_of[fn+i];
    fh->twin() = Halffacet_of[fn+(i/2*2)+((i+1)%2)];
    fh->boundary_entry_objects().push_back(make_object(SEdge_of[sen+bnd[i]]));
    fh->incident_volume() = Volume_of[((i%4) == 1 || (i%4 == 2)) ? 1 : 0];
    if(i<4) {
      hz = i % 2 ? -1 : 1;
      hx = hy = 0;
    }
    else if(i<8) {
      hy = i % 2 ? -1 : 1;
      hx = hz = 0;
    }
    else {
      hx = i % 2 ? -1 : 1;
      hz = hy = 0;
    }
    hw = ((i%4) == 1 || (i%4) == 2) ? 1 : -1;
    fh->plane() = Infi_box::create_extended_plane(hx,hy,hz,hw);
    fh->mark() = 1;
  }

  Volume_of[0]->shell_entry_objects().push_back(make_object(SFace_of[sfn]));
  Volume_of[0]->mark() = 0;
  Volume_of[1]->shell_entry_objects().push_front(make_object(SFace_of[sfn+1]));

  int sprevOff[6] = {4,3,0,5,2,1};
  int snextOff[6] = {2,5,4,1,0,3};
  int prevIdx[48] = {7,12,15,26,29,10,
                     1,18,21,32,35,4,
                     19,0,3,38,41,22,
                     13,6,9,44,47,16,
                     31,36,39,2,5,34,
                     25,42,45,8,11,28,
                     43,24,27,14,17,46,
                     37,30,33,20,23,40};
  int nextIdx[48] = {13,6,27,14,11,28,
                     19,0,33,20,5,34,
                     1,18,39,2,23,40,
                     7,12,45,8,17,46,
                     37,30,3,38,35,4,
                     43,24,9,44,29,10,
                     25,42,15,26,47,16,
                     31,36,21,32,41,22};
  int factIdx[48] = {1,0,9,8,5,4,
                     0,1,11,10,4,5,
                     0,1,8,9,7,6,
                     1,0,10,11,6,7,
                     3,2,8,9,4,5,
                     2,3,10,11,5,4,
                     2,3,9,8,6,7,
                     3,2,11,10,7,6};
  int sgn[24] = {1,1,1,-1,1,-1,
                 -1,-1,1,1,-1,-1,
                 1,-1,-1,-1,-1,1,
                 -1,1,-1,1,1,1};

  for(int i = 0; i < 48; ++i) {
    SHalfedge_handle seh = SEdge_of[sen+i];

    seh->twin() = SEdge_of[sen+(i/2*2)+((i+1)%2)];
    seh->sprev() = SEdge_of[sen+sprevOff[i%6]+(i/6*6)];
    seh->snext() = SEdge_of[sen+snextOff[i%6]+(i/6*6)];
    seh->source() = Edge_of[en+((i+1)%6)/2+(i/6)*3];
    seh->incident_sface() = SFace_of[sfn+(i%2)+(i/6)*2];
    seh->prev() = SEdge_of[sen+prevIdx[i]];
    seh->next() = SEdge_of[sen+nextIdx[i]];
    seh->facet() = Halffacet_of[fn+factIdx[i]];
    if(i%6 < 2) {
      hz = (i%2) ? sgn[i/2] * (-1) : sgn[i/2];
      hx = hy = 0;
    }
    else if(i%6 < 4) {
      hx = (i%2) ? sgn[i/2] * (-1) : sgn[i/2];
      hz = hy = 0;
    }
    else {
      hy = (i%2) ? sgn[i/2] * (-1) : sgn[i/2];
      hx = hz = 0;
    }
    seh->circle() = Sphere_circle(Plane_3(RT(hx),RT(hy),RT(hz),RT(0)));
    seh->mark() = 1;
  }

  int volIdx[8] = {0,1,1,0,1,0,0,1};

  for(int i = 0; i < 16; ++i) {
    SFace_handle sfh = SFace_of[sfn+i];
    sfh->center_vertex() = Vertex_of[vn+(i/2)];
    sfh->boundary_entry_objects().push_back(make_object(SEdge_of[sen+(i/2*6)+(i%2)]));
    this->sncp()->store_sm_boundary_item(SEdge_of[sen+(i/2*6)+(i%2)],
                                          --(sfh->sface_cycles_end()));
    int cIdx = i%2 ? 1-volIdx[i/2] : volIdx[i/2];
    sfh->volume() = Volume_of[cIdx];
    sfh->mark() = cIdx ? Volume_of[1]->mark() : 0;
  }
}

} //namespace CGAL
#endif //CGAL_SNC_IO_PARSER_H
