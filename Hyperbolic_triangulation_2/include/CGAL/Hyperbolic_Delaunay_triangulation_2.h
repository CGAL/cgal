// Copyright (c) 2010-2018  INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mikhail Bogdanov
//                 Monique Teillaud <Monique.Teillaud@inria.fr>

#ifndef CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
#define CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/license/Hyperbolic_triangulation_2.h>

#include <CGAL/Hyperbolic_triangulation_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <algorithm>
#include <stack>
#include <set>
#include <vector>

namespace CGAL {

template < class Gt,
           class Tds = Triangulation_data_structure_2 <
                         Triangulation_vertex_base_2<Gt>,
                         Hyperbolic_triangulation_face_base_2<Gt> > >
class Hyperbolic_Delaunay_triangulation_2
  : private Delaunay_triangulation_2<Gt,Tds>
{
public:
  typedef Hyperbolic_Delaunay_triangulation_2<Gt, Tds>          Self;
  typedef Delaunay_triangulation_2<Gt,Tds>                      Base;

  typedef typename Tds::size_type                               size_type;
  typedef typename Tds::Vertex_handle                           Vertex_handle;
  typedef typename Tds::Face_handle                             Face_handle;
  typedef typename Tds::Edge                                    Edge;

  typedef Gt                                                    Geom_traits;
  typedef typename Geom_traits::FT                              FT;
  typedef typename Geom_traits::Hyperbolic_point_2              Point;
  typedef typename Geom_traits::Hyperbolic_Voronoi_point_2      Hyperbolic_Voronoi_point;
  typedef typename Geom_traits::Hyperbolic_segment_2            Hyperbolic_segment;
  typedef typename Geom_traits::Hyperbolic_triangle_2           Hyperbolic_triangle;

  // Tag to distinguish regular triangulations from others
  typedef Tag_false                                             Weighted_tag;

  // Tag to distinguish periodic triangulations from others
  typedef Tag_false                                             Periodic_tag;

#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Base::cw;
  using Base::ccw;
  using Base::geom_traits;
#endif

  enum Locate_type {
    VERTEX = 0,
    EDGE,
    FACE,
    OUTSIDE_CONVEX_HULL,
    OUTSIDE_AFFINE_HULL
  };

  typedef typename Geom_traits::Side_of_oriented_hyperbolic_segment_2 Side_of_oriented_hyperbolic_segment;
  typedef typename Geom_traits::Is_Delaunay_hyperbolic                Is_Delaunay_hyperbolic;

  Hyperbolic_Delaunay_triangulation_2(const Geom_traits& gt = Geom_traits())
    : Delaunay_triangulation_2<Gt,Tds>(gt), _gt(gt)
  {
  }

  Hyperbolic_Delaunay_triangulation_2(const Hyperbolic_Delaunay_triangulation_2<Gt,Tds> &tr)
    : Delaunay_triangulation_2<Gt,Tds>(tr), _gt()
  {
    CGAL_postcondition(this->is_valid());
  }

  template<class InputIterator>
  Hyperbolic_Delaunay_triangulation_2(InputIterator first, InputIterator last,
                                      const Geom_traits& gt = Geom_traits())
    : Delaunay_triangulation_2<Gt,Tds>(gt), _gt(gt)
  {
    insert(first, last);
  }

  /*************************************
      Circulators and iterators
  *************************************/
private:
  // This class is used to generate the iterators.
  class Non_hyperbolic_tester
  {
    const Self *t;

  public:
    Non_hyperbolic_tester() {} // needs a default constructor for Filter_iterator
    Non_hyperbolic_tester(const Self *tr) : t(tr) {}

    bool operator()(const typename Base::All_vertices_iterator & vit) const { return Base::is_infinite(vit); }
    bool operator()(const typename Base::All_faces_iterator & fit) const { return !t->is_Delaunay_hyperbolic(fit); }
    bool operator()(const typename Base::All_edges_iterator & eit) const
    {
      return !t->is_Delaunay_hyperbolic(eit->first, eit->second);
    }
  };

  Non_hyperbolic_tester non_hyperbolic_tester() const { return Non_hyperbolic_tester(this); }

public:
  class Hyperbolic_faces_iterator
    : public Filter_iterator<typename Base::All_faces_iterator, Non_hyperbolic_tester>
  {
    typedef Filter_iterator<typename Base::All_faces_iterator, Non_hyperbolic_tester> PBase;
    typedef Hyperbolic_faces_iterator                                                 Self;

  public:
    Hyperbolic_faces_iterator() : PBase() {}
    Hyperbolic_faces_iterator(const PBase &b) : PBase(b) {}
    Self & operator++() { PBase::operator++(); return *this; }
    Self & operator--() { PBase::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
    operator const Face_handle() const { return PBase::base(); }
  };

  Hyperbolic_faces_iterator hyperbolic_faces_begin() const
  {
    if(dimension() < 2)
      return hyperbolic_faces_end();

    return CGAL::filter_iterator(Base::all_faces_end(),
                                 Non_hyperbolic_tester(this),
                                 Base::all_faces_begin());
  }

  Hyperbolic_faces_iterator hyperbolic_faces_end() const
  {
    return CGAL::filter_iterator(Base::all_faces_end(), Non_hyperbolic_tester(this));
  }

  typedef Filter_iterator<typename Base::All_edges_iterator, Non_hyperbolic_tester> Hyperbolic_edges_iterator;

  Hyperbolic_edges_iterator hyperbolic_edges_begin() const
  {
    if(dimension() < 1)
      return hyperbolic_edges_end();

    return CGAL::filter_iterator(Base::all_edges_end(),
                                 Non_hyperbolic_tester(this),
                                 Base::all_edges_begin());
  }

  Hyperbolic_edges_iterator hyperbolic_edges_end() const
  {
    return CGAL::filter_iterator(Base::all_edges_end(), Non_hyperbolic_tester(this));
  }


  template <typename HTriangulation>
  class Hyperbolic_adjacent_vertex_circulator
    : Base::Vertex_circulator
  {
    typedef typename Base::Vertex_circulator        VBase;
    typedef Hyperbolic_adjacent_vertex_circulator   Self;
    typedef typename Tds::Vertex                    Vertex;

    Vertex_handle _v;
    Face_handle   pos;
    int _ri;
    int _iv;

  public:
    Hyperbolic_adjacent_vertex_circulator(const HTriangulation& tri) : VBase(), _v(Vertex_handle()), pos(Face_handle()), _ri(0), _tri(tri) {}

    Hyperbolic_adjacent_vertex_circulator(Vertex_handle v, const HTriangulation& tri, Face_handle fh = Face_handle())
      : VBase(v, fh), _tri(tri)
    {
      _v = v;
      if (fh == Face_handle())
        pos = _v->face();
      else
        pos = fh;
      _iv = pos->index(_v);

      bool ok = false;
      do
      {
        _ri = cw(_iv);
        if (!_tri.is_Delaunay_hyperbolic(pos, ccw(_iv)))
        {
          _ri = ccw(_iv);
          if (!_tri.is_Delaunay_hyperbolic(pos, cw(_iv)))
          {
            pos = pos->neighbor(cw(_iv));
            _iv = pos->index(_v);
          }
          else
          {
            ok = true;
          }
        }
        else {
          ok = true;
        }
      } while (!ok);

    }

    Self& operator++()
    {
      pos = pos->neighbor(cw(_iv));
      _iv = pos->index(_v);

      bool ok = false;
      do
      {
        _ri = cw(_iv);
        if (!_tri.is_Delaunay_hyperbolic(pos, ccw(_iv)))
        {
          _ri = ccw(_iv);
          if (!_tri.is_Delaunay_hyperbolic(pos, cw(_iv)))
          {
            pos = pos->neighbor(cw(_iv));
            _iv = pos->index(_v);
          }
          else
          {
            ok = true;
          }
        }
        else {
          ok = true;
        }
      } while (!ok);

      return *this;
    }


    Self& operator--()
    {
      pos = pos->neighbor(ccw(_iv));
      _iv = pos->index(_v);

      bool ok = false;
      do
      {
        _ri = ccw(_iv);
        if (!_tri.is_Delaunay_hyperbolic(pos, cw(_iv)))
        {
          _ri = cw(_iv);
          if (!_tri.is_Delaunay_hyperbolic(pos, ccw(_iv)))
          {
            pos = pos->neighbor(ccw(_iv));
            _iv = pos->index(_v);
          }
          else
          {
            ok = true;
          }
        }
        else {
          ok = true;
        }
      } while (!ok);

      return *this;
    }

    bool operator==(const Self &vc) const
    {
      return (this->_v == vc->_v &&
              this->pos == vc->pos &&
              this->_ri == vc->_ri &&
              this->_iv == vc->_iv);
    }

    bool operator!=(const Self &vc) const { return !this->operator==(vc); }
    bool operator==(const Vertex_handle &vh) const { return (this->pos->vertex(_ri) == vh); }
    bool operator!=(const Vertex_handle &vh) const { return !this->operator==(vh); }

    bool is_empty() const { return (this->pos == Face_handle() && this->_v == Vertex_handle()); }

    Vertex& operator*() const
    {
      CGAL_precondition(pos != Face_handle() && _v != Vertex_handle());
      return *(pos->vertex(_ri));
    }

    Vertex* operator->() const
    {
      CGAL_precondition(pos != Face_handle() && _v != Vertex_handle());
      return &*(pos->vertex(_ri));
    }

    Vertex_handle base() const { return pos->vertex(_ri); }
    operator Vertex_handle() const { return pos->vertex(_ri); }

  private:
    const HTriangulation& _tri;
  };

private:
  Geom_traits _gt;

public:
  Tds& tds()
  {
    return Base::tds();
  }

  const Tds& tds() const
  {
    return Base::tds();
  }

  Geom_traits& geom_traits()
  {
    return _gt;
  }

  const Geom_traits& geom_traits() const
  {
    return _gt;
  }

  void swap (Self&  tr)
  {
    Base::swap(tr);

    Geom_traits t = _gt;
    _gt = tr._gt;
    tr._gt = t;

    this->mark_finite_non_hyperbolic_faces();
    tr.mark_finite_non_hyperbolic_faces();
  }

  Self& operator=(const Self &tr)
  {
    Self newone = Self(tr);
    this->swap(newone);

    return *this;
  }

  bool operator==(const Self& tr )
  {
    if (tr.number_of_vertices() != this->number_of_vertices())
      return false;

    if (tr.number_of_hyperbolic_faces() != this->number_of_hyperbolic_faces())
      return false;

    if (tr.number_of_vertices() == 0)
      return true; // as in Periodic_2_triangulation_2

    return Base::operator==(tr);

  }


  typedef Hyperbolic_adjacent_vertex_circulator<Self> Vertex_circulator;
  typedef typename Tds::Edge_circulator         Edge_circulator;
  typedef typename Tds::Face_circulator         Face_circulator;

  typedef Hyperbolic_faces_iterator             All_faces_iterator;
  All_faces_iterator all_faces_begin() const { return hyperbolic_faces_begin(); }
  All_faces_iterator all_faces_end() const { return hyperbolic_faces_end(); }

  typedef Hyperbolic_edges_iterator             All_edges_iterator;
  All_edges_iterator all_edges_begin() const { return hyperbolic_edges_begin(); }
  All_edges_iterator all_edges_end() const { return hyperbolic_edges_end(); }

  typedef typename Base::Finite_vertices_iterator All_vertices_iterator;
  All_vertices_iterator all_vertices_begin() const { return Base::finite_vertices_begin(); }
  All_vertices_iterator all_vertices_end() const { return Base::finite_vertices_end(); }

  // The declarations below are required by apply_to_range: do not document!
  typedef All_vertices_iterator Finite_vertices_iterator;
  Finite_vertices_iterator finite_vertices_begin() const { return all_vertices_begin(); }
  Finite_vertices_iterator finite_vertices_end() const { return all_vertices_end(); }
  typedef All_faces_iterator Finite_faces_iterator;
  Finite_faces_iterator finite_faces_begin() const { return all_faces_begin(); }
  Finite_faces_iterator finite_faces_end() const { return all_faces_end(); }

  // Algebraic_kernel_for_circles_2 needs this for some reason
  typedef typename Base::Line_face_circulator     Line_face_circulator;

  void clear() { Base::clear(); }

  template<class OutputItFaces>
  OutputItFaces find_conflicts(const Point& p, OutputItFaces fit, Face_handle start = Face_handle()) const
  {
    return Base::get_conflicts(p, fit, start);
  }

  Vertex_handle insert(const Point& p,
                       Face_handle start = Face_handle())
  {
    Vertex_handle v = Base::insert(p, start);
    mark_star_faces(v);
    ensure_hyperbolic_face_handle(v);

    return v;
  }

  // Note that this is _Base::Locate_type_ and not Locate_type.
  // Do _not_ use this function with the result of a call to HDT2::locate(), which is a locate
  // on the hyperbolic triangulation, and not the underlying Delaunay triangulation (in which
  // new points are inserted, and thus it needs to be DT2::locate_type).
  Vertex_handle insert(const Point& p,
                       typename Base::Locate_type lt,
                       Face_handle loc, int li)
  {
    Vertex_handle v = Base::insert(p, lt, loc, li);
    mark_star_faces(v);
    ensure_hyperbolic_face_handle(v);

    return v;
  }

#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  template < class InputIterator >
  std::ptrdiff_t insert(InputIterator first, InputIterator last,
                        std::enable_if_t<
                          std::is_base_of<Point, typename std::iterator_traits<InputIterator>::value_type>::value
                        >* = nullptr)
#else
  template < class InputIterator >
  std::ptrdiff_t insert(InputIterator first, InputIterator last)
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  {
    size_type n = Base::insert(first, last);

    mark_finite_non_hyperbolic_faces();
    for(All_vertices_iterator vit = all_vertices_begin(); vit != all_vertices_end(); ++vit)
      ensure_hyperbolic_face_handle(vit);

    return n;
  }

  void remove(Vertex_handle v)
  {
    CGAL_precondition(tds().is_vertex(v));
    std::vector<Vertex_handle> nbr;
    bool dim_was_2 = false;
    if (this->dimension() == 2)
    {
      dim_was_2 = true;
      typename Base::Vertex_circulator nbv = Base::incident_vertices(v), done(nbv);
      do
      {
        nbr.push_back(nbv);
      } while (++nbv != done);
    }

    Base::remove(v);

    if (dim_was_2)
    {
      for (unsigned int i = 0; i < nbr.size(); ++i)
      {
        mark_star_faces(nbr[i]);
        ensure_hyperbolic_face_handle(nbr[i]);
      }
    }
  }

  template <class VertexRemoveIterator>
  void remove(VertexRemoveIterator first, VertexRemoveIterator last)
  {
    for (VertexRemoveIterator vit = first; vit != last; ++vit)
    {
      remove(*vit);
    }
  }

  template <typename T>
  bool is_infinite(T v) const { return Base::is_infinite(v); }
  bool is_infinite(Face_handle f, int i) const { return Base::is_infinite(f, i); }

  bool is_Delaunay_hyperbolic(Face_handle f) const
  {
    if(dimension() <= 1)
      return false;

    return f->hyperbolic_data().is_Delaunay_hyperbolic();
  }

  bool is_Delaunay_hyperbolic(Face_handle f, int i) const
  {
    if(dimension() <= 1)
      return false;

    if(is_infinite(f, i))
      return false;

    if(f->hyperbolic_data().is_Delaunay_non_hyperbolic(i))
      return false;

    // another incident face and corresponding index
    Face_handle f2 = f->neighbor(i);
    int i2 = f2->index(f);

    if(f2->hyperbolic_data().is_Delaunay_non_hyperbolic(i2))
      return false;

    return true;
  }

  bool is_Delaunay_hyperbolic(const Edge& e) const
  {
    return is_Delaunay_hyperbolic(e.first, e.second);
  }

  bool is_Delaunay_hyperbolic(const Edge_circulator& ec) const
  {
    return is_Delaunay_hyperbolic(*ec);
  }

  bool is_Delaunay_hyperbolic(const All_edges_iterator& ei) const
  {
    return is_Delaunay_hyperbolic(*ei);
  }

private:
  /*
    During the insertion of a new point in the triangulation, the added vertex points to a face.
    This function ensures that the face to which the vertex points is hyperbolic (if there exists one).
  */
  void ensure_hyperbolic_face_handle(Vertex_handle v)
  {
    if(dimension() <= 1)
      return;

    Face_circulator fc = this->incident_faces(v), done(fc);
    if(fc != 0)
    {
      do
      {
        if(is_Delaunay_hyperbolic(fc))
        {
          v->set_face(fc);
          break;
        }
      }
      while(++fc != done);
    }
  }

  Oriented_side side_of_hyperbolic_triangle(const Point& p, const Point& q, const Point& r,
                                            const Point& query, Locate_type &lt, int& li) const
  {

    // The triangle (p,q,r) must be Delaunay hyperbolic
    CGAL_precondition(geom_traits().is_Delaunay_hyperbolic_2_object()(p, q, r));
    CGAL_precondition(query != p && query != q && query != r);

    // Point p is assumed to be at index 0, q at index 1 and r at index 2 in the face.
    li = -1;

    if(query == p)
    {
      lt = VERTEX;
      li = 0;
      return ON_ORIENTED_BOUNDARY;
    }

    if(query == q)
    {
      lt = VERTEX;
      li = 1;
      return ON_ORIENTED_BOUNDARY;
    }

    if(query == r)
    {
      lt = VERTEX;
      li = 2;
      return ON_ORIENTED_BOUNDARY;
    }

    Oriented_side cp1 = geom_traits().side_of_oriented_hyperbolic_segment_2_object()(p, q, query);
    if(cp1 == ON_ORIENTED_BOUNDARY)
    {
      lt = EDGE;
      li = 2;
      return ON_ORIENTED_BOUNDARY;
    }

    Oriented_side cp2 = geom_traits().side_of_oriented_hyperbolic_segment_2_object()(q, r, query);
    if(cp2 == ON_ORIENTED_BOUNDARY)
    {
      lt = EDGE;
      li = 0;
      return ON_ORIENTED_BOUNDARY;
    }

    Oriented_side cp3 = geom_traits().side_of_oriented_hyperbolic_segment_2_object()(r, p, query);
    if(cp3 == ON_ORIENTED_BOUNDARY)
    {
      lt = EDGE;
      li = 1;
      return ON_ORIENTED_BOUNDARY;
    }

    Oriented_side cs1 = geom_traits().side_of_oriented_hyperbolic_segment_2_object()(p, q, r);
    Oriented_side cs2 = geom_traits().side_of_oriented_hyperbolic_segment_2_object()(q, r, p);
    Oriented_side cs3 = geom_traits().side_of_oriented_hyperbolic_segment_2_object()(r, p, q);

    // Cannot be on the boundary here.
    lt = FACE;
    li = 4;

    if(cs1 != cp1 || cs2 != cp2 || cs3 != cp3)
      return ON_NEGATIVE_SIDE;
    else
      return ON_POSITIVE_SIDE;
  }

  void mark_finite_non_hyperbolic_faces() const
  {
    if(dimension() <= 1)
      return;

    for(auto fit = Base::all_faces_begin(); fit != Base::all_faces_end(); ++fit)
      fit->hyperbolic_data().set_Delaunay_hyperbolic(); // finite & hyperbolic

    Face_handle ifh = Base::infinite_face();
    ifh->hyperbolic_data().set_infinite();

    std::stack<Face_handle> to_visit;
    to_visit.push(ifh);

    std::set<Face_handle> visited_faces; // @todo squat tds_data()

    while(!to_visit.empty())
    {
      Face_handle fh = to_visit.top();
      to_visit.pop();

      if(!visited_faces.insert(fh).second) // already visited previously
        continue;

      for(int i = 0; i<3; ++i)
      {
        Face_handle nfh = fh->neighbor(i);
        mark_face(nfh);

        if(is_Delaunay_hyperbolic(nfh))
          continue;

        to_visit.push(nfh);
      }
    }
  }

  void mark_star_faces(Vertex_handle v) const
  {
    if(dimension() <= 1)
      return;

    Face_handle f = v->face();
    Face_handle start(f);
    do
    {
      mark_face(f);

      int i = f->index(v);
      f = f->neighbor(ccw(i));
    }
    while(f != start);
  }

  void mark_face(const Face_handle f) const
  {
    if(is_infinite(f))
    {
      f->hyperbolic_data().set_infinite();
    }
    else
    {
      int idx;
      bool flag = geom_traits().is_Delaunay_hyperbolic_2_object()(point(f,0),
                                                                  point(f,1),
                                                                  point(f,2),
                                                                  idx);

      if(flag)
        f->hyperbolic_data().set_Delaunay_hyperbolic(); // finite & hyperbolic
      else
        f->hyperbolic_data().set_Delaunay_non_hyperbolic(idx); // finite but not hyperbolic
    }
  }

public:
  Line_face_circulator line_walk(const Point& p, const Point& q, Face_handle f = Face_handle()) const
  {
    return Base::line_walk(p, q, f);
  }

  Hyperbolic_triangle hyperbolic_triangle(const Face_handle f) const
  {
    CGAL_precondition(!is_infinite(f));
    return Base::triangle(f);
  }

  // needed by DT_2: do not document!
  Hyperbolic_triangle triangle(const Face_handle f) const
  {
    CGAL_precondition(!is_infinite(f));
    return hyperbolic_triangle(f);
  }

  Hyperbolic_segment hyperbolic_segment(const Face_handle f, const int i) const
  {
    return geom_traits().construct_hyperbolic_segment_2_object()(point(f,cw(i)),
                                                                 point(f,ccw(i)));
  }

  Hyperbolic_segment hyperbolic_segment(const Edge& e) const
  {
    Face_handle f = e.first;
    int i = e.second;
    return hyperbolic_segment(f, i);
  }

  Hyperbolic_segment hyperbolic_segment(const Edge_circulator& e) const { return hyperbolic_segment(*e); }

  // needed by DT_2: do not document!
  Hyperbolic_segment segment(const Face_handle f, const int i) const { return hyperbolic_segment(f,i); }
  Hyperbolic_segment segment(const Edge& e) const { return hyperbolic_segment(e); }
  Hyperbolic_segment segment(const Edge_circulator& e) const { return hyperbolic_segment(e); }

  size_type number_of_vertices() const { return Base::number_of_vertices(); }
  Vertex_circulator adjacent_vertices(Vertex_handle v) const { return Vertex_circulator(v, *this); }

  size_type number_of_hyperbolic_faces() const
  {
    return std::distance(hyperbolic_faces_begin(), hyperbolic_faces_end());
  }

  size_type number_of_hyperbolic_edges() const
  {
    return std::distance(hyperbolic_edges_begin(), hyperbolic_edges_end());
  }

  int dimension() const { return Base::dimension(); }

  Hyperbolic_Voronoi_point dual(Face_handle f) const
  {
    CGAL_precondition(is_Delaunay_hyperbolic(f));
    return geom_traits().construct_hyperbolic_circumcenter_2_object()(point(f,0),
                                                                      point(f,1),
                                                                      point(f,2));
  }

  Hyperbolic_segment dual(const Edge& e) const { return dual(e.first, e.second); }

  Hyperbolic_segment dual(Face_handle f, int i) const
  {
    CGAL_precondition(is_Delaunay_hyperbolic(f, i));

    if(dimension() == 1)
    {
      const Point& p = point(f,cw(i));
      const Point& q = point(f,ccw(i));

      // hyperbolic line
      Hyperbolic_segment line = geom_traits().construct_hyperbolic_bisector_2_object()(p, q);
      return line;
    }

    Face_handle n = f->neighbor(i);
    int in = n->index(f);

    bool fhyp = is_Delaunay_hyperbolic(f);
    bool nhyp = is_Delaunay_hyperbolic(n);

    // both faces are non_hyperbolic, but the common edge is hyperbolic
    if(!fhyp && !nhyp)
    {
      const Point& p = point(f,ccw(i));
      const Point& q = point(f,cw(i));

      // hyperbolic line
      Hyperbolic_segment line = geom_traits().construct_hyperbolic_bisector_2_object()(p, q);
      return line;
    }

    // both faces are hyperbolic
    if(fhyp && nhyp)
    {
      const Point& p = point(f,ccw(i));
      const Point& q = point(f,cw(i));

      Hyperbolic_segment s = geom_traits().construct_hyperbolic_bisector_2_object()(
                               p, q, point(f,i), point(n,in));
      return s;
    }

    // one of the incident faces is non_hyperbolic
    Face_handle hyp_face = f;

    if(!fhyp)
    {
      hyp_face = n;
      i = in;
    }

    const Point& p = point(hyp_face,ccw(i));
    const Point& q = point(hyp_face,cw(i));

    Hyperbolic_segment ray = geom_traits().construct_hyperbolic_bisector_2_object()(
                               p, q, point(hyp_face,i));
    return ray;
  }

public:
  const Point& point(const Vertex_handle vh) const
  {
    CGAL_precondition(!is_infinite(vh));
    return vh->point();
  }

  const Point& point(const Face_handle fh, const int i) const
  {
    CGAL_precondition(!is_infinite(fh->vertex(i)));
    CGAL_precondition(0 <= i && i <= 2);
    return fh->vertex(i)->point();
  }

  Point& point(const Vertex_handle vh)
  {
    CGAL_precondition(!is_infinite(vh));
    return vh->point();
  }

  Point& point(const Face_handle fh, const int i)
  {
    CGAL_precondition(!is_infinite(fh->vertex(i)));
    CGAL_precondition(0 <= i && i <= 2);
    return fh->vertex(i)->point();
  }

  bool is_valid()
  {
    if (!Base::is_valid())
      return false;

    for (Hyperbolic_faces_iterator fit = hyperbolic_faces_begin(); fit != hyperbolic_faces_end(); fit++)
    {
      if (!is_Delaunay_hyperbolic(fit))
        return false;
    }

    for (Hyperbolic_edges_iterator eit = hyperbolic_edges_begin(); eit != hyperbolic_edges_end(); eit++)
    {
      if (!is_Delaunay_hyperbolic(eit))
        return false;
    }

    return true;
  }

  Face_handle locate(const Point& p, const Face_handle hint = Face_handle()) const
  {
    Locate_type lt;
    int li;
    return locate(p, lt, li, hint);
  }

  Face_handle locate(const Point& query, Locate_type& lt, int &li, Face_handle hint = Face_handle()) const
  {
    // Perform a Euclidean location first and get close to the hyperbolic face containing the query point
    typename Base::Locate_type blt;
    Face_handle fh = Base::locate(query, blt, li, hint);

    switch(blt)
    {
      case Base::VERTEX: lt = VERTEX; break;
      case Base::EDGE: lt = EDGE; break;
      case Base::FACE: lt = FACE; break;
      case Base::OUTSIDE_CONVEX_HULL: lt = OUTSIDE_CONVEX_HULL; break;
      case Base::OUTSIDE_AFFINE_HULL: lt = OUTSIDE_AFFINE_HULL; break;
    }

    if(lt == VERTEX)
      return fh;

    if(lt == OUTSIDE_CONVEX_HULL || lt == OUTSIDE_AFFINE_HULL)
      return Face_handle();

    CGAL_assertion(!is_infinite(fh));

    // This case corresponds to when the point is located on a Euclidean edge.
    if(lt == EDGE)
    {
      // Here because the call to `side_of_hyperbolic_triangle` might change `li`
      Face_handle mfh = fh->neighbor(li);

      if(is_Delaunay_hyperbolic(fh))
      {
        const Point& p = point(fh, 0);
        const Point& q = point(fh, 1);
        const Point& r = point(fh, 2);

        Oriented_side side = side_of_hyperbolic_triangle(p, q, r, query, lt, li);
        if(side != ON_NEGATIVE_SIDE)
          return fh;
      }

      if(is_Delaunay_hyperbolic(mfh))
      {
        const Point& p = point(mfh, 0);
        const Point& q = point(mfh, 1);
        const Point& r = point(mfh, 2);
        Oriented_side side = side_of_hyperbolic_triangle(p, q, r, query, lt, li);

        if(side != ON_NEGATIVE_SIDE) {
          return fh;
        } else {
          lt = OUTSIDE_CONVEX_HULL;
          return Face_handle();
        }
      }
      else
      {
        lt = OUTSIDE_CONVEX_HULL;
        return Face_handle();
      }
    }

    // Here, the face has been located in the Euclidean face fh
    const Point& p = point(fh, 0);
    const Point& q = point(fh, 1);
    const Point& r = point(fh, 2);

    if(!is_Delaunay_hyperbolic(fh))
    {
      // Need to check if the point lies on one of the sides of the face
      // Note that at least one side is Delaunay hyperbolic!
      if(geom_traits().side_of_oriented_hyperbolic_segment_2_object()(p,q,query) == ON_ORIENTED_BOUNDARY)
      {
        lt = EDGE;
        li = 2;
        return fh;
      }
      else if(geom_traits().side_of_oriented_hyperbolic_segment_2_object()(q,r,query) == ON_ORIENTED_BOUNDARY)
      {
        lt = EDGE;
        li = 0;
        return fh;
      }
      else if(geom_traits().side_of_oriented_hyperbolic_segment_2_object()(r,p,query) == ON_ORIENTED_BOUNDARY)
      {
        lt = EDGE;
        li = 1;
        return fh;
      }

      lt = OUTSIDE_CONVEX_HULL;
      return Face_handle();
    }

    Oriented_side side = side_of_hyperbolic_triangle(p, q, r, query, lt, li);
    if(side != ON_NEGATIVE_SIDE) {
      return fh;
    } else {
      // Here, the point lies in a face that is a neighbor to fh
      for(int i = 0; i < 3; ++i) {
        Face_handle nfh = fh->neighbor(i);
        if(is_Delaunay_hyperbolic(nfh))
        {
          Oriented_side nside = side_of_hyperbolic_triangle(point(nfh,0),
                                                            point(nfh,1),
                                                            point(nfh,2),
                                                            query, lt, li);
          if(nside != ON_NEGATIVE_SIDE)
            return nfh;
        }
      }

      // At this point, the point lies outside of the convex hull of the triangulation,
      // since it has not been found in any of the hyperbolic faces adjacent to fh.
      lt = OUTSIDE_CONVEX_HULL;
      return Face_handle();
    }

    // We never reach this point, but we have to make the compiler happy
    lt = OUTSIDE_CONVEX_HULL;
    return Face_handle();
  }
};


} // end namespace CGAL

#endif // CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
