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

#include <stack>
#include <set>

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
    CGAL_triangulation_postcondition(this->is_valid());
  }

  template<class InputIterator>
  Hyperbolic_Delaunay_triangulation_2(InputIterator first, InputIterator last,
                                      const Geom_traits& gt = Geom_traits())
    : Delaunay_triangulation_2<Gt,Tds>(gt), _gt(gt)
  {
    insert(first, last);
    for(All_vertices_iterator vit=all_vertices_begin(); vit!=all_vertices_end(); ++vit)
      ensure_hyperbolic_face_handle(vit);
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
      Edge e(eit->first, eit->second);
      return !t->is_Delaunay_hyperbolic(e);
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
        if (_tri.is_finite_non_hyperbolic(pos, ccw(_iv)))
        {
          _ri = ccw(_iv);
          if (_tri.is_finite_non_hyperbolic(pos, cw(_iv)))
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
        if (_tri.is_finite_non_hyperbolic(pos, ccw(_iv)))
        {
          _ri = ccw(_iv);
          if (_tri.is_finite_non_hyperbolic(pos, cw(_iv)))
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
        if (_tri.is_finite_non_hyperbolic(pos, cw(_iv)))
        {
          _ri = cw(_iv);
          if (_tri.is_finite_non_hyperbolic(pos, ccw(_iv)))
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
      CGAL_triangulation_precondition(pos != Face_handle() && _v != Vertex_handle());
      return *(pos->vertex(_ri));
    }

    Vertex* operator->() const
    {
      CGAL_triangulation_precondition(pos != Face_handle() && _v != Vertex_handle());
      return &*(pos->vertex(_ri));
    }

    Vertex_handle base() const { return pos->vertex(_ri); }
    operator Vertex_handle() const { return pos->vertex(_ri); }

  private:
    const HTriangulation& _tri;
  };

private:
  Geom_traits   _gt;

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

  void mark_star(Vertex_handle v) const
  {
    if(!is_star_bounded(v))
      mark_star_faces(v);
  }

  template<class OutputItFaces>
  OutputItFaces find_conflicts(const Point& p, OutputItFaces fit, Face_handle start = Face_handle()) const
  {
    return Base::get_conflicts(p, fit, start);
  }

  Vertex_handle insert(const Point& p,
                       Face_handle start = Face_handle())
  {
    Vertex_handle v = Base::insert(p, start);
    mark_star(v);
    ensure_hyperbolic_face_handle(v);

    return v;
  }

  Vertex_handle insert(const Point& p,
                       typename Base::Locate_type lt,
                       Face_handle loc, int li)
  {
    Vertex_handle v = Base::insert(p, lt, loc, li);
    mark_star(v);
    ensure_hyperbolic_face_handle(v);

    return v;
  }

#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  template < class InputIterator >
  std::ptrdiff_t insert(InputIterator first, InputIterator last,
                        typename boost::enable_if<
                          boost::is_base_of<Point, typename std::iterator_traits<InputIterator>::value_type>
                        >::type* = nullptr)
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
    CGAL_triangulation_precondition(tds().is_vertex(v));
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



  /*
    Needed by DT_2: do not document!
  */
  template <typename T>
  bool is_infinite(T v) const { return Base::is_infinite(v); }

  bool is_Delaunay_hyperbolic(Face_handle f) const
  {
    return !Base::is_infinite(f) && !is_finite_non_hyperbolic(f);
  }

  bool is_Delaunay_hyperbolic(Face_handle f, int i) const
  {
    return !Base::is_infinite(f, i) && !is_finite_non_hyperbolic(f, i);
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
  class Face_data
  {
  private:
    // a finite face is non_hyperbolic if its circumscribing circle intersects the circle at infinity
    bool _is_Delaunay_hyperbolic;

    // defined only if the face is finite and non_hyperbolic
    unsigned int _non_hyperbolic_edge;

  public:
    Face_data() : _is_Delaunay_hyperbolic(true), _non_hyperbolic_edge(UCHAR_MAX) {}

    unsigned int get_non_hyperbolic_edge() const
    {
      CGAL_triangulation_precondition(!_is_Delaunay_hyperbolic);
      CGAL_triangulation_precondition(_non_hyperbolic_edge <= 2);

      return _non_hyperbolic_edge;
    }

    void set_non_hyperbolic_edge(unsigned int uschar)
    {
      CGAL_triangulation_precondition(!_is_Delaunay_hyperbolic);
      CGAL_triangulation_precondition(uschar <= 2);

      _non_hyperbolic_edge = uschar;
    }

    bool get_is_Delaunay_hyperbolic() const { return _is_Delaunay_hyperbolic; }
    void set_is_Delaunay_hyperbolic(bool flag) { _is_Delaunay_hyperbolic = flag; }
  };

  /*
    During the insertion of a new point in the triangulation, the added vertex points to a face.
    This function ensures that the face to which the vertex points is hyperbolic.
  */
  void ensure_hyperbolic_face_handle(Vertex_handle v)
  {
    if(dimension() > 2)
    {
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
      CGAL_triangulation_postcondition(is_Delaunay_hyperbolic(v->face()));
    }
  }

  Oriented_side side_of_hyperbolic_triangle(const Point& p, const Point& q, const Point& r,
                                            const Point& query, Locate_type &lt, int& li) const
  {
    // The triangle (p,q,r) must be Delaunay hyperbolic
    CGAL_triangulation_precondition(geom_traits().is_Delaunay_hyperbolic_2_object()(p, q, r));

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
    if(cs1 != cp1 || cs2 != cp2 || cs3 != cp3)
      return ON_NEGATIVE_SIDE;
    else
      return ON_POSITIVE_SIDE;
  }

  int get_finite_non_hyperbolic_edge(Face_handle f) const
  {
    CGAL_triangulation_precondition(is_finite_non_hyperbolic(f));
    Face_data fd = object_cast<Face_data>(f->tds_data());
    return fd.get_non_hyperbolic_edge();
  }

  bool is_finite_non_hyperbolic(Face_handle f) const
  {
    if(const Face_data* td = object_cast<Face_data>(&f->tds_data()))
    {
      return !td->get_is_Delaunay_hyperbolic();
    }
    else
    {
      return false;
    }
  }

  bool is_finite_non_hyperbolic(Face_handle f, int i) const
  {
    if(dimension() <= 1)
      return false;

    if(is_finite_non_hyperbolic(f) && get_finite_non_hyperbolic_edge(f) == i)
      return true;

    // another incident face and corresponding index
    Face_handle f2 = f->neighbor(i);
    int i2 = f2->index(f);

    if(is_finite_non_hyperbolic(f2) && get_finite_non_hyperbolic_edge(f2) == i2)
      return true;

    return false;
  }

  bool is_finite_non_hyperbolic(const Edge& e) const
  {
    return is_finite_non_hyperbolic(e.first, e.second);
  }

  // Depth-first search (dfs) and marking the finite non_hyperbolic faces.
  void mark_finite_non_hyperbolic_faces() const
  {
    if(dimension() <= 1)
      return;

    std::set<Face_handle> visited_faces;

    // maintain a stack to be able to backtrack
    // to the most recent faces which neighbors are not visited
    std::stack<Face_handle> backtrack;

    // start from a face with infinite vertex
    Face_handle current = Base::infinite_face();

    // mark it as visited
    visited_faces.insert(current);

    // put the element whose neighbors we are going to explore.
    backtrack.push(current);

    Face_handle next;

    while(!backtrack.empty())
    {
      // take a face
      current = backtrack.top();

      // start visiting the neighbors
      int i = 0;
      for(; i<3; ++i)
      {
        next = current->neighbor(i);

        // if a neighbor is already visited, then stop going deeper
        if(visited_faces.find(next) != visited_faces.end())
          continue;

        visited_faces.insert(next);
        mark_face(next);

        // go deeper if the neighbor is non_hyperbolic
        if(!is_Delaunay_hyperbolic(next))
        {
          backtrack.push(next);
          break;
        }
      }

      // if all the neighbors are already visited, then remove "current" face.
      if(i == 3)
        backtrack.pop();
    }
  }

  // check if a star is bounded by finite faces
  bool is_star_bounded(Vertex_handle v) const
  {
    if(dimension() <= 1)
      return true;

    Face_handle f = v->face();
    Face_handle next;
    int i;
    Face_handle start(f);
    Face_handle opposite_face;

    do
    {
      i = f->index(v);
      next = f->neighbor(ccw(i));  // turn ccw around v

      opposite_face = f->neighbor(i);
      if(!is_Delaunay_hyperbolic(opposite_face))
        return false;

      f = next;
    }
    while(next != start);

    return true;
  }


  void mark_star_faces(Vertex_handle v) const
  {
    if(dimension() <= 1)
      return;

    Face_handle f = v->face();
    Face_handle start(f), next;
    int i;
    do
    {
      i = f->index(v);
      next = f->neighbor(ccw(i));  // turn ccw around v

      mark_face(f);

      f = next;
    } while(next != start);
  }

  void mark_face(const Face_handle f) const
  {
    Is_Delaunay_hyperbolic del;
    int idx;
    bool flag = del(point(f,0),
                    point(f,1),
                    point(f,2),
                    idx);

    Face_data fd;
    fd.set_is_Delaunay_hyperbolic(flag);

    if(!flag)
      fd.set_non_hyperbolic_edge(idx);

    f->tds_data() = make_object(fd);
  }

public:

  Line_face_circulator line_walk(const Point& p, const Point& q, Face_handle f = Face_handle()) const
  {
    return Base::line_walk(p, q, f);
  }

  Hyperbolic_triangle hyperbolic_triangle(const Face_handle f) const { return Base::triangle(f); }

  // needed by DT_2: do not document!
  Hyperbolic_triangle triangle(const Face_handle f) const { return hyperbolic_triangle(f); }

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
    CGAL_triangulation_precondition(is_Delaunay_hyperbolic(f));
    return geom_traits().construct_hyperbolic_circumcenter_2_object()(point(f,0),
                                                                      point(f,1),
                                                                      point(f,2));
  }

  Hyperbolic_segment dual(const Edge& e) const { return dual(e.first, e.second); }

  Hyperbolic_segment dual(Face_handle f, int i) const
  {
    CGAL_triangulation_precondition(is_Delaunay_hyperbolic(f, i));

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

    // both faces are non_hyperbolic, but the incident edge is hyperbolic
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

  const Point point(const Vertex_handle vh) const
  {
    return vh->point();
  }

  const Point point(const Face_handle fh, const int i) const
  {
    CGAL_triangulation_precondition(0 <= i);
    CGAL_triangulation_precondition(i <= 2);
    return fh->vertex(i)->point();
  }


  Point point(const Vertex_handle vh)
  {
    return vh->point();
  }

  Point point(const Face_handle fh, const int i)
  {
    CGAL_triangulation_precondition(0 <= i);
    CGAL_triangulation_precondition(i <= 2);
    return fh->vertex(i)->point();
  }



  bool is_valid()
  {
    if (Base::is_valid())
    {
      for (Hyperbolic_faces_iterator fit = hyperbolic_faces_begin(); fit != hyperbolic_faces_end(); fit++)
      {
        if (!is_Delaunay_hyperbolic(fit))
        {
          return false;
        }
      }
      for (Hyperbolic_edges_iterator eit = hyperbolic_edges_begin(); eit != hyperbolic_edges_end(); eit++)
      {
        if (!is_Delaunay_hyperbolic(eit))
        {
          return false;
        }
      }
      return true;
    }
    return false;
  }

  Face_handle locate(const Point& p, const Face_handle hint = Face_handle()) const
  {
    Locate_type lt;
    int li;
    return locate(p, lt, li, hint);
  }

  Face_handle locate(const Point& query, Locate_type& lt, int &li, Face_handle hint = Face_handle()) const
  {
    // Perform an Euclidean location first and get close to the hyperbolic face containing the query point
    typename Base::Locate_type blt;
    Face_handle fh = Base::locate(query, blt, li, hint);

    if(blt == Base::VERTEX) {
      lt = VERTEX;
    } else {
      if(blt == Base::EDGE) {
        lt = EDGE;
      } else {
        if(blt == Base::FACE) {
          lt = FACE;
        } else {
          if(blt == Base::OUTSIDE_CONVEX_HULL) {
            lt = OUTSIDE_CONVEX_HULL;
          } else {
            lt = OUTSIDE_AFFINE_HULL;
          }
        }
      }
    }

    if(lt == VERTEX)
      return fh;

    if(lt == OUTSIDE_CONVEX_HULL || lt == OUTSIDE_AFFINE_HULL)
      return Face_handle();

    // This case corresponds to when the point is located on an Euclidean edge.
    if(lt == EDGE)
    {
      Point p = point(fh, 0);
      Point q = point(fh, 1);
      Point r = point(fh, 2);

      if(geom_traits().is_Delaunay_hyperbolic_2_object()(p, q, r))
      {
        Oriented_side side = side_of_hyperbolic_triangle(p, q, r, query, lt, li);
        if(side == ON_ORIENTED_BOUNDARY) {
          lt = EDGE;
          return fh;
        } else {
          if(side == ON_POSITIVE_SIDE) {
            lt = FACE;
            return fh;
          } else {
            // do nothing -- we still have to check the neighboring face
          }
        }
      }

      p = point(fh, ccw(li));
      q = point(Base::mirror_vertex(fh, li));
      r = point(fh, cw(li));

      if(geom_traits().is_Delaunay_hyperbolic_2_object()(p, q, r))
      {
        Oriented_side side = side_of_hyperbolic_triangle(p, q, r, query, lt, li);
        if(side == ON_ORIENTED_BOUNDARY) {
          lt = EDGE;
          return fh;
        } else {
          if(side == ON_POSITIVE_SIDE) {
            lt = FACE;
            return fh;
          } else {
            // There is nothing to be done now -- the point is outside the convex hull of the triangulation
            lt = OUTSIDE_CONVEX_HULL;
            return Face_handle();
          }
        }
      }
    }

    // Here, the face has been located in the Euclidean face lh
    const Point& p = point(fh, 0);
    const Point& q = point(fh, 1);
    const Point& r = point(fh, 2);
    int idx;
    if(!geom_traits().is_Delaunay_hyperbolic_2_object()(p, q, r, idx))
    {
      // Need to check if the point lies on one of the sides of the face
      // Note that at least one side is Delaunay hyperbolic!
      if (geom_traits().side_of_oriented_hyperbolic_segment_2_object()(p,q,query) == ON_ORIENTED_BOUNDARY ||
          geom_traits().side_of_oriented_hyperbolic_segment_2_object()(q,r,query) == ON_ORIENTED_BOUNDARY ||
          geom_traits().side_of_oriented_hyperbolic_segment_2_object()(r,p,query) == ON_ORIENTED_BOUNDARY   )
          lt = EDGE;
      else
        lt = OUTSIDE_CONVEX_HULL;
      return Face_handle();
    }

    Oriented_side side = side_of_hyperbolic_triangle(p, q, r, query, lt, li);
    if(side == ON_POSITIVE_SIDE) {
      lt = FACE;
      return fh;
    } else {
      if(side == ON_ORIENTED_BOUNDARY) {
        lt = EDGE;
        return fh;
      } else {
        // Here, the point lies in a face that is a neighbor to fh
        for(int i = 0; i < 3; i++) {
          Face_handle nfh = fh->neighbor(i);
          if(geom_traits().is_Delaunay_hyperbolic_2_object()(point(nfh,0),
                                                             point(nfh,1),
                                                             point(nfh,2)))
          {
            Oriented_side nside = side_of_hyperbolic_triangle(point(nfh,0),
                                                              point(nfh,1),
                                                              point(nfh,2),
                                                              query, lt, li);
            if(nside == ON_POSITIVE_SIDE) {
              lt = FACE;
              return nfh;
            } else if(nside == ON_ORIENTED_BOUNDARY) {
              lt = EDGE;
              return nfh;
            }
          }
        }

        // At this point, the point lies outside of the convex hull of the triangulation,
        // since it has not been found in any of the hyperbolic faces adjacent to fh.
        lt = OUTSIDE_CONVEX_HULL;
        return Face_handle();
      }
    }

    // We never reach this point, but we have to make the compiler happy
    lt = OUTSIDE_CONVEX_HULL;
    return Face_handle();
  }
};


} // end namespace CGAL

#endif // CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
