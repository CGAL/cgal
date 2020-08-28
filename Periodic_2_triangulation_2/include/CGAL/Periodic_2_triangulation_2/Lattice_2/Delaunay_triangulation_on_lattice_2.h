// Copyright(c) 2020 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Mael Rouxel-Labbé
//                 Georg Osang

// See also publication:
//   Georg Osang, Mael Rouxel-Labbé, and Monique Teillaud
//   Generalizing CGAL Periodic Delaunay Triangulations
//   ESA 2020
//   https://doi.org/10.4230/LIPIcs.ESA.2020.75

#ifndef CGAL_P2T2_DELAUNAY_TRIANGULATION_ON_LATTICE_2_H
#define CGAL_P2T2_DELAUNAY_TRIANGULATION_ON_LATTICE_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>

#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Triangulation_face_base_on_lattice_2.h>
#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Triangulation_iterators_on_lattice_2.h>
#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Triangulation_circulators_on_lattice_2.h>
#include <CGAL/Periodic_2_triangulation_2/Lattice_2/Triangulation_on_lattice_2.h>
#include <CGAL/Periodic_2_triangulation_2/Square_flat_torus_2/Delaunay_triangulation_on_square_flat_torus_2.h>
#include <CGAL/Periodic_2_triangulation_vertex_base_2.h>
#include <CGAL/Periodic_2_triangulation_2/IO/periodic_2_triangulation_2_io.h>
#include <CGAL/draw_triangulation_2.h>

#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/utility.h>

#include <utility>
#include <iostream>
#include <fstream>
#include <string>

namespace CGAL {

template <typename Gt_,
          typename Tds_ = Triangulation_data_structure_2<
                            Periodic_2_triangulation_vertex_base_2<Gt_>,
                            Triangulation_face_base_on_lattice_2<Gt_> > >
class Delaunay_triangulation_on_lattice_2
  : public Triangulation_on_lattice_2<
             Gt_, Tds_,
             Delaunay_triangulation_2<Gt_, Tds_>,
             Delaunay_triangulation_on_square_flat_torus_2<Gt_, Tds_> >
{
public:
  typedef Gt_                                                             Geom_traits;
  typedef Tds_                                                            Triangulation_data_structure;

private:
  typedef Geom_traits                                                     Gt;
  typedef Triangulation_data_structure                                    Tds;
  typedef Delaunay_triangulation_2<Gt, Tds>                               DT2;
  typedef Delaunay_triangulation_on_square_flat_torus_2<Gt, Tds>          P2DT2;

  typedef Triangulation_on_lattice_2<Gt, Tds, DT2, P2DT2>                 Base;
  typedef Delaunay_triangulation_on_lattice_2<Gt, Tds>                    Self;

public:
  typedef typename Base::Offset                                           Offset;
  typedef typename Base::Domain                                           Domain;
  typedef typename Base::Covering_sheets                                  Covering_sheets;

public:
  typedef typename Base::Vertex                                           Vertex;
  typedef typename Base::Vertex_handle                                    Vertex_handle;
  typedef typename Base::Edge                                             Edge;
  typedef typename Base::Face_handle                                      Face_handle;
  typedef typename Base::Edge_circulator                                  Edge_circulator;
  typedef typename Base::Edge_iterator                                    Edge_iterator;

  typedef typename Base::FT                                               FT;
  typedef typename Base::Point                                            Point;
  typedef typename Base::Segment                                          Segment;

  typedef typename Base::Periodic_point                                   Periodic_point;
  typedef typename Base::Locate_type                                      Locate_type;

public:
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2
  using Base::cw;
  using Base::ccw;

  using Base::dimension;
  using Base::domain;
  using Base::geom_traits;
  using Base::tds;
  using Base::t2;
  using Base::p2t2;

  using Base::canonical_vertex;
  using Base::construct_segment;
  using Base::get_neighbor_offset;
  using Base::get_canonical_face;
  using Base::is_canonical;
  using Base::is_1_cover;
  using Base::neighbor;
  using Base::periodic_point;
#endif

  // This is only used when inserting the canonical point in the multiple cover setting (DT2)
  class Conflict_tester
  {
    const Self& _tr_ptr;
    Point _p;

  public:
    /// Constructor
    Conflict_tester(const Point& pt, const Self& tr_ptr) : _tr_ptr(tr_ptr), _p(pt) { }

    template <typename OutputIterator>
    void operator()(OutputIterator out) const
    {
      _tr_ptr.t2().get_conflicts(_p, out);
    }

    void set_point(const Point& pt) { _p = pt; }
    const Point& point() const { return _p; }
  };

  class Cover_manager
  {
    Self& tr;

  public:
    Cover_manager(Self& tr) : tr(tr) { }

    void insert_unsatisfying_element(const Face_handle fh) const {
      tr.insert_face_with_too_big_circumradius(fh);
    }

    void delete_unsatisfying_element(const Face_handle fh) const {
      tr.delete_face_with_too_big_circumradius(fh);
    }

    bool can_be_converted_to_1_sheet() const {
      return tr.can_be_converted_to_1_sheet();
    }

    // needed for Tr::remove()
//    bool update_cover_data_during_management(Face_handle new_f,
//                                             const std::vector<Face_handle>& new_faces,
//                                             const bool abort_if_cover_change)
//    {
//      return tr.update_cover_data_during_management(new_f, new_faces, abort_if_cover_change);
//    }
  };

private:
  FT sq_circumradius_threshold;
  std::set<Face_handle> faces_with_too_big_circumradius;

public:
  /// Constructors
  Delaunay_triangulation_on_lattice_2(const Geom_traits& gt)
    :  Base(gt)
  {
    update_cover_data_after_setting_domain(gt.get_domain());
  }

  Delaunay_triangulation_on_lattice_2(const Domain& lattice)
    : Delaunay_triangulation_on_lattice_2(Geom_traits(lattice))
  { }

  template <class InputIterator>
  Delaunay_triangulation_on_lattice_2(InputIterator first, InputIterator beyond,
                                      const Geom_traits& gt)
    : Delaunay_triangulation_on_lattice_2(gt)
  {
    insert(first, beyond);
  }

  template <class InputIterator>
  Delaunay_triangulation_on_lattice_2(InputIterator first, InputIterator beyond,
                                      const Domain& lattice)
    : Delaunay_triangulation_on_lattice_2(first, beyond, Geom_traits(lattice))
  { }

  void swap(Delaunay_triangulation_on_lattice_2& tr)
  {
    std::swap(sq_circumradius_threshold, tr.sq_circumradius_threshold);
    std::swap(faces_with_too_big_circumradius, tr.faces_with_too_big_circumradius);
    Base::swap(static_cast<Base&>(tr));
  }

  // @todo (virtual vertices cannot be naively copied)
  Delaunay_triangulation_on_lattice_2& operator=(const Delaunay_triangulation_on_lattice_2&) = delete;
  Delaunay_triangulation_on_lattice_2(const Delaunay_triangulation_on_lattice_2&) = delete;
  void copy_triangulation(const Delaunay_triangulation_on_lattice_2 &tr) = delete;

  void clear()
  {
    Base::clear();
    clear_covering_data();
    faces_with_too_big_circumradius.clear();
  }

  /// Domain setting
public:
  void update_cover_data_after_setting_domain(const Domain& domain)
  {
    sq_circumradius_threshold = 0.0625 * domain.systole_sq_length();
  }

  void clear_covering_data()
  {
    faces_with_too_big_circumradius.clear();
  }

  void set_domain(const Domain& domain)
  {
    clear();
    Base::set_domain(domain);
    update_cover_data_after_setting_domain(domain());
  }

  /// Access
public:
  size_t too_big_faces() const { return faces_with_too_big_circumradius.size(); }

  /// Predicates
public:
  Oriented_side side_of_oriented_circle(const Point& p0, const Point& p1, const Point& p2, const Point& p,
                                        bool perturb) const
  {
    if(is_1_cover())
      return p2t2().side_of_oriented_circle(p0, p1, p2, p, perturb);
    else
      return t2().side_of_oriented_circle(p0, p1, p2, p, perturb);
  }

  Oriented_side side_of_oriented_circle(const Face_handle f, const Point& p,
                                        bool perturb = false) const
  {
    if(is_1_cover())
      return p2t2().side_of_oriented_circle(f, p, perturb);
    else
      return t2().side_of_oriented_circle(f, p, perturb);
  }

  FT compare_squared_circumradius_to_threshold(const Face_handle f, const FT threshold)
  {
    if(is_1_cover())
    {
      const Periodic_point p0 = periodic_point(f, 0);
      const Periodic_point p1 = periodic_point(f, 1);
      const Periodic_point p2 = periodic_point(f, 2);

      return geom_traits().compare_squared_radius_2_object()(p0.first, p1.first, p2.first,
                                                             p0.second, p1.second, p2.second,
                                                             threshold);
    }
    else
    {
      return geom_traits().compare_squared_radius_2_object()(t2().point(f, 0),
                                                             t2().point(f, 1),
                                                             t2().point(f, 2),
                                                             threshold);
    }
  }

public:
  Vertex_handle nearest_vertex(const Point& p,
                               Face_handle f = Face_handle())
  {
    const Point cp = geom_traits().get_domain().construct_canonical_point(p);

    if(is_1_cover())
    {
      return p2t2().nearest_vertex(cp, f);
    }
    else
    {
      Vertex_handle nv = t2().nearest_vertex(cp, f);
      return canonical_vertex(nv);
    }
  }

  /// Conflicts
public:
  template<class OutputItFaces , class OutputItBoundaryEdges >
  std::pair<OutputItFaces, OutputItBoundaryEdges>
  get_conflicts_and_boundary(const Point& p,
                             OutputItFaces fit,
                             OutputItBoundaryEdges eit,
                             Face_handle start) const
  {
    CGAL_assertion(false); // @todo
  }

  template<class OutputItFaces >
  OutputItFaces get_conflicts(const Point& p,
                              OutputItFaces fit,
                              Face_handle start) const
  {
    CGAL_assertion(false); // @todo
    return fit;
  }

  template<class OutputItBoundaryEdges >
  OutputItBoundaryEdges get_boundary_of_conflicts(const Point &p,
                                                  OutputItBoundaryEdges eit,
                                                  Face_handle start) const
  {
    CGAL_assertion(false); // @todo
    return eit;
  }

public:
  /// Insertion and removal
  Vertex_handle insert(const Point& p,
                       Locate_type lt,
                       Face_handle f,
                       int li)
  {
    CGAL_triangulation_precondition(is_canonical(p));

    Conflict_tester conflict_tester(p, *this);
    Cover_manager cover_manager(*this);

    return Base::insert_in_conflict(p, lt, f, li, conflict_tester, cover_manager);
  }

  Vertex_handle insert(const Point& p)
  {
    const Point cp = geom_traits().get_domain().construct_canonical_point(p);

    Conflict_tester conflict_tester(p, *this);
    Cover_manager cover_manager(*this);

    return Base::insert_in_conflict(cp, conflict_tester, cover_manager);
  }

  Vertex_handle push_back(const Point& p)
  {
    return insert(p);
  }

  template <class InputIterator>
  void insert(InputIterator first, InputIterator beyond)
  {
    // @todo spatial sorting
    while(first != beyond)
      insert(*first++);

    CGAL_triangulation_postcondition(is_valid());
  }

  void remove(Vertex_handle)
  {
    // @todo
    CGAL_triangulation_assertion(false);
  }

  Vertex_handle move_if_no_collision(Vertex_handle v, const Point& p)
  {
    Locate_type lt;
    int li;
    Vertex_handle inserted;
    Face_handle loc = locate(p, lt, li, v->face());

    if(lt == Base::VERTEX)
      return v;
    else
      // This can be optimized by checking whether we can move v->point() to p
      return insert(p, lt, loc, li);
  }

  Vertex_handle move_point(Vertex_handle v, const Point& p)
  {
    if(v->point() == p)
      return v;

    Vertex_handle w = move_if_no_collision(v, p);
    if(w != v)
    {
      remove(v);
      return w;
    }
    return v;
  }

private:
  bool can_be_converted_to_1_sheet() const
  {
    return faces_with_too_big_circumradius.empty();
  }

  void insert_face_with_too_big_circumradius(const Face_handle fh)
  {
    if(compare_squared_circumradius_to_threshold(fh, sq_circumradius_threshold) != CGAL::SMALLER)
      faces_with_too_big_circumradius.insert(fh);
  }

  void delete_face_with_too_big_circumradius(const Face_handle fh)
  {
    faces_with_too_big_circumradius.erase(fh);
  }

  /// Duals
public:
  Point circumcenter(Face_handle f) const
  {
    CGAL_triangulation_precondition(is_canonical(f));

    if(is_1_cover())
      return p2t2().circumcenter(f);
    else
      return t2().circumcenter(f);
  }

  Point dual(const Face_handle f) const
  {
    CGAL_triangulation_precondition(is_canonical(f));

    if(is_1_cover())
      return p2t2().dual(f);
    else
      return t2().dual(f);
  }

  Segment dual(const Edge& e) const
  {
    CGAL_triangulation_precondition(is_canonical(e));

    if(is_1_cover())
    {
      return p2t2().dual(e);
    }
    else
    {
      Face_handle f = e.first;
      int i = e.second;

      if(is_canonical(f))
      {
        Face_handle nb = neighbor(f, i);
        CGAL_assertion(is_canonical(nb));
        Point p0 = dual(f);
        Point p1 = dual(nb);
        Offset o = get_neighbor_offset(f, i);
        return construct_segment(p0, p1, o, Offset());
      }
      else
      {
        // Find the edge again such that e.first is a canonical face
        Vertex_handle v = canonical_vertex(f->vertex(ccw(i)));
        f = get_canonical_face(f);
        for(int j=0; j<3; ++j)
        {
          if(v == canonical_vertex(f->vertex(j)))
          {
            i = cw(j);
            break;
          }
        }

        return dual(Edge(f, i));
      }
    }
  }

  Segment dual(const Edge_circulator& ec) const { return dual(*ec); }
  Segment dual(const Edge_iterator& ei) const { return dual(*ei); }

  template<class Stream>
  Stream& draw_dual(Stream& ps)
  {
    Edge_iterator eit = Base::edges_begin(), eend = Base::edges_end();
    for(; eit!=eend; ++eit)
    {
      Segment d = dual(eit);
      ps << "2 " << d.source() << " 0 " << d.target() << " 0\n";
    }

    return ps;
  }

  /// is_valid and IO
public:
  bool is_valid(const bool verbose = false, int level = 0)
  {
    // Note that checking the base here also checks Delaunay validity
    // because the triangulations t2 and p2t2 are Delaunay triangulations
    bool result = Base::is_valid(verbose, level);

    if(is_1_cover())
    {
      if(!faces_with_too_big_circumradius.empty())
      {
        result = false;
        if(verbose)
        {
          std::cerr << "Error: periodic triangulation in one cover "
                    << "but there are faces with too big circumradius" << std::endl;
        }

        return false;
      }
    }

    return true;
  }

  template <typename GT, typename TDS>
  friend std::ostream& operator<<(std::ostream& os, const Delaunay_triangulation_on_lattice_2<GT, TDS>& tr)
  {
    if(tr.is_1_cover())
      CGAL::write_PD2T2_to_OFF(os, tr.p2t2());
    else
      CGAL::write_DT2_to_OFF(os, tr.t2());

    return os;
  }
};

} // namespace CGAL

#endif // CGAL_P2T2_DELAUNAY_TRIANGULATION_ON_LATTICE_2_H
