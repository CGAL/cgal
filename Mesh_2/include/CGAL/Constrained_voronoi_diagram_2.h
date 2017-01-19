// Copyright (c) 2013 INRIA Sophia-Antipolis (France),
//               2014-2015 GeometryFactory (France).
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
// Author(s) : Jane Tournois, Raul Gallegos, Pierre Alliez
//

#ifndef CGAL_CONSTRAINED_VORONOI_DIAGRAM_2_H
#define CGAL_CONSTRAINED_VORONOI_DIAGRAM_2_H

#include <CGAL/license/Mesh_2.h>


#include <utility>
#include <stack>
#include <CGAL/iterator.h>
#include <CGAL/tuple.h>
#include <CGAL/Kernel/global_functions_2.h>

namespace CGAL {

template <class Cdt>
class Cvd_cell_2
{
  typedef typename Cdt::Vertex_handle    Vertex_handle;

public:
  typedef typename Cdt::Geom_traits::Segment_2    Segment;
  typedef typename Cdt::Geom_traits::Ray_2        Ray;
  typedef typename Cdt::Geom_traits::Point_2      Point;
  typedef CGAL::Dispatch_output_iterator<
    CGAL::cpp11::tuple<Segment, Ray>,
    CGAL::cpp11::tuple<std::back_insert_iterator<std::vector<Segment> >,
                       std::back_insert_iterator<std::vector<Ray> > >
    > Construction_dispatcher;

private:
  Vertex_handle m_vertex; //generator
  std::vector<Segment> m_segments;
  std::vector<Ray>     m_rays;
  bool m_is_valid;
  
public:
  Cvd_cell_2(Vertex_handle v)
    : m_vertex(v)
    , m_segments()
    , m_rays()
    , m_is_valid(false)
  { }

  bool operator<(const Cvd_cell_2& cell) const
  {
    return m_vertex < cell.vertex();
  };

  Vertex_handle vertex() const { return m_vertex; }

  bool is_valid() const { return m_is_valid; }
  bool& is_valid()      { return m_is_valid; }

  bool is_infinite() const 
  {
    return !m_rays.empty();
  }
  bool is_empty() const
  {
    return m_rays.empty() && m_segments.empty();
  }

public:
  //construction iterators
  std::back_insert_iterator<std::vector<Segment> > segment_output_iterator()
  {
    return std::back_inserter(m_segments);
  }
  std::back_insert_iterator<std::vector<Ray> > ray_output_iterator()
  {
    return std::back_inserter(m_rays);
  }

public:
  std::size_t number_of_vertices() const
  {
    return m_segments.size();
  }

  Point point(const std::size_t& i) const
  {
    CGAL_assertion(i >= 0 && i < m_segments.size());
    return m_segments[i].source();
  }

public:
  //access iterators
  typedef typename std::vector<Segment>::iterator  segment_iterator;
  typedef typename std::vector<Ray>::iterator      ray_iterator;
  typedef typename std::vector<Segment>::const_iterator  const_segment_iterator;
  typedef typename std::vector<Ray>::const_iterator      const_ray_iterator;

  segment_iterator segments_begin()        { return m_segments.begin(); }
  segment_iterator segments_end()          { return m_segments.end(); }
  const_segment_iterator segments_cbegin() const { return m_segments.cbegin();}
  const_segment_iterator segments_cend()   const { return m_segments.cend(); }

  ray_iterator rays_begin()       { return m_rays.begin(); }
  ray_iterator rays_end()         { return m_rays.end(); }
  const_ray_iterator rays_cbegin() const { return m_rays.cbegin();}
  const_ray_iterator rays_cend()   const { return m_rays.cend(); }

  CGAL_assertion_code(
public:
  bool is_simply_ccw_oriented() const
  {
    typedef typename Cdt::Geom_traits::Vector_2 Vector;
    const_segment_iterator sit = segments_cbegin();
    Segment s1 = *sit++;
    for(; sit != segments_cend(); ++sit)
    {
      Segment s2 = *sit;
      if(s1.target() != s2.source())
        return false;

      Point p = vertex()->point();
      Vector v1(p, s1.source());
      Vector v2(p, s1.target());
      Vector v3(p, s2.target());

      if(CGAL::orientation(v1, v2) != CGAL::LEFT_TURN
        || CGAL::orientation(v2, v3) != CGAL::LEFT_TURN)
        return false;

      s1 = s2;
    }
    return true;
  }
  );//end CGAL_assertion_code

}; //end CLASS Cvd_cell_2


// Cdt should be of the type Constrained_Delaunay_triangulation_2
// and the face base shoul be Constrained_Delaunay_triangulation_face_base_2
template <class Cdt>
class Constrained_voronoi_diagram_2
{
public:
  typedef Constrained_voronoi_diagram_2<Cdt>         Cvd;
  typedef Cvd_cell_2<Cdt>                            Cvd_cell;
  typedef typename Cvd_cell::Construction_dispatcher Construction_dispatcher;

public:
  // typedefs for basic primitives 
  typedef typename Cdt::Geom_traits       Geom_traits;
  typedef typename Cdt::Intersection_tag Intersection_tag;

  typedef typename Cdt::Constraint             Constraint;
  typedef typename Cdt::Vertex_handle          Vertex_handle;
  typedef typename Cdt::Face_handle            Face_handle;
  typedef typename Cdt::Edge                   Edge;
  typedef typename Cdt::All_faces_iterator     All_faces_iterator;
  typedef typename Cdt::Finite_faces_iterator  Finite_faces_iterator;
  typedef typename Cdt::Finite_edges_iterator  Finite_edges_iterator;
  typedef typename Cdt::Finite_vertices_iterator Finite_vertices_iterator;
  typedef typename Cdt::Face_circulator        Face_circulator;
  typedef typename Cdt::Edge_circulator        Edge_circulator;
  typedef typename Cdt::size_type              size_type;
  typedef typename Cdt::Locate_type            Locate_type;

  typedef typename Geom_traits::FT                   FT;
  typedef typename Geom_traits::Point_2              Point;
  typedef typename Geom_traits::Vector_2             Vector;
  typedef typename Geom_traits::Line_2               Line;
  typedef typename Cdt::Segment                 Segment;
  typedef typename Cdt::Triangle                Triangle;

protected:
  const Cdt& m_cdt;

public:
  Constrained_voronoi_diagram_2(const Cdt& cdt)
    : m_cdt(cdt)
  {
  }

  //----------------------------------------------------------------
  //--------------------ABOUT FACES SIGHT---------------------------
  //----------------------------------------------------------------

public:
  // blind = false IFF each face sees its circumcenter
  void tag_all_faces_blind(const bool blind) 
  {
    for(All_faces_iterator f = m_cdt.all_faces_begin();
         f != m_cdt.all_faces_end();
         ++f)
      f->set_blind(blind);
  }

  // blind test for each face
  // if true, set corresponding barrier constraint
  void tag_faces_blind()
  {
    if(m_cdt.dimension() < 2)
      return;

    tag_all_faces_blind(false);

    // for each constrained edge, mark blinded triangles
    for(Finite_edges_iterator e = m_cdt.finite_edges_begin();
         e != m_cdt.finite_edges_end();
         ++e)
    {
      Edge edge = *e;
      if(m_cdt.is_constrained(edge))
      {
        tag_neighbors_blind(edge);
        tag_neighbors_blind(m_cdt.mirror_edge(edge));
      }
    }
  }

private:
  // test face for blindness with respect to the edge constraint
  void tag_face_blind(Face_handle& f, const Edge& constraint)
  {  
    if(segment_hides_circumcenter(m_cdt.segment(constraint),
                                  m_cdt.triangle(f)))
    {
      f->set_blind(true);
      f->set_blinding_constraint(constraint);
    }
  }

  // predicate: returns true if the triangle tr and its circumcenter
  // are on the opposite side of the segment seg
  bool segment_hides_circumcenter(const Segment& seg,
                                  const Triangle& tr)
  {
    Point a = seg.source();
    Point b = seg.target();
    double dX = b.x() - a.x();
    double dY = b.y() - a.y();

    const Point& p0 = tr[0];
    const Point& p1 = tr[1];
    const Point& p2 = tr[2];
    double R0 = p0.x()*p0.x() + p0.y()*p0.y();
    double R1 = p1.x()*p1.x() + p1.y()*p1.y();
    double R2 = p2.x()*p2.x() + p2.y()*p2.y();
    double denominator = (p1.x()-p0.x())*(p2.y()-p0.y()) +
                                   (p0.x()-p2.x())*(p1.y()-p0.y());

    double det = 2*denominator * (a.x()*dY - a.y()*dX)
                    - (R2-R1) * (p0.x()*dX + p0.y()*dY)
                    - (R0-R2) * (p1.x()*dX + p1.y()*dY)
                    - (R1-R0) * (p2.x()*dX + p2.y()*dY);
    return (det <= 0);
  }

  // tags with their sights, with respect to the Edge constraint,
  // seed and its neighbor faces, on the same side of Edge than seed.
  void tag_neighbors_blind(const Edge& constraint)
  {
    CGAL_assertion(m_cdt.is_constrained(constraint));
    Face_handle seed = constraint.first;

    if(!m_cdt.is_infinite(seed) 
       && !seed->is_blind() 
       && !m_cdt.triangle(seed).is_degenerate() )
       //to avoid flat triangles outside the domain
    {
      std::stack<Face_handle> faces;
      faces.push(seed);

      while(!faces.empty())
      {
        Face_handle f = faces.top();
        faces.pop();
        this->tag_face_blind(f, constraint);
        if(f->is_blind())
          this->push_unvisited_neighbors(f, faces);
      }
    }
  }

  // puts in the stack the unvisited (un-tagged) neighbor faces of f
  void push_unvisited_neighbors(const Face_handle& f,
                                std::stack<Face_handle>& faces) const
  {
    for(int i=0; i<3; ++i)
    {
      Face_handle fi = f->neighbor(i);
      Edge edge_i = Edge(f, i);
      if(!m_cdt.is_constrained(edge_i) &&
          !fi->is_blind() &&
          !m_cdt.is_infinite(fi)) 
        faces.push(fi);
    }
  }

  /*--------------------------------------------------------------
  ---------------------- BVD CONSTRUCTION ------------------------
  --------------------------------------------------------------*/

public:
  Cvd_cell cvd_cell(Vertex_handle v) const
  {
    Cvd_cell cell(v);
    
    typename Cvd_cell::Construction_dispatcher oit =
      CGAL::dispatch_output<typename Cvd_cell::Segment,
                            typename Cvd_cell::Ray>(
                              cell.segment_output_iterator(),
                              cell.ray_output_iterator());
    cvd_cell(v, oit);
    cell.is_valid() = true;

    return cell;
  }

// assemble a cell of the bounded Voronoi diagram
// incident to vertex v
// OutputIterator should be able to collect Segments and Rays
private:
  template<typename OutputIterator>
  OutputIterator cvd_cell(Vertex_handle v, OutputIterator oit) const
  {
    if(bvd_cell_is_infinite(v))
      infinite_cvd_cell(v, oit);
    else
      finite_cvd_cell(v, oit);
    return oit;
  }

private:
  template <typename OutputIterator>
  OutputIterator finite_cvd_cell(Vertex_handle v, OutputIterator oit) const
  {
    std::vector<Point> polygon;
    
    CGAL_assertion(!m_cdt.is_infinite(v));
    Face_circulator face = m_cdt.incident_faces(v);
    Face_circulator end = face;
    Face_circulator next = face;

    CGAL_For_all(face, end)
    {
      next++;
      Line line(m_cdt.circumcenter(face), m_cdt.circumcenter(next));
      Point intersect;

      if(!face->is_blind()) //face sees
      {
        polygon.push_back(m_cdt.circumcenter(face));
        if(next->is_blind())  //next doesn't
        {
          CGAL_assertion(do_intersect(line, m_cdt.segment(next->blinding_constraint())));
          CGAL::assign(intersect,
            CGAL::intersection(line, Line(m_cdt.segment(next->blinding_constraint()))));
          polygon.push_back(intersect);
        }
      }
      else //face doesn't see
      {
        if(!next->is_blind()) //next sees
        {
          CGAL_assertion(do_intersect(line, m_cdt.segment(face->blinding_constraint())));
          CGAL::assign(intersect,
            CGAL::intersection(line, Line(m_cdt.segment(face->blinding_constraint()))));
          polygon.push_back(intersect);
        }
        else //next doesn't
        {
          if(face->blinding_constraint() != next->blinding_constraint()
            && face->blinding_constraint() != m_cdt.mirror_edge(next->blinding_constraint()))
            // the 2 blinding_constraints are different
          {
            CGAL_assertion(do_intersect(line, m_cdt.segment(face->blinding_constraint())));
            CGAL::assign(intersect,
              CGAL::intersection(line, Line(m_cdt.segment(face->blinding_constraint()))));
            polygon.push_back(intersect);

            Point intersection2;
            CGAL_assertion(do_intersect(line, m_cdt.segment(next->blinding_constraint())));
            CGAL::assign(intersection2,
              CGAL::intersection(line, Line(m_cdt.segment(next->blinding_constraint()))));
            polygon.push_back(intersection2);
          }
          //else: it's the same constraint--> do nothing
        }
      }
    }//end CGAL_For_all

    std::size_t nbp = polygon.size();
    for(std::size_t i = 0; i < nbp; ++i)
      *oit++ = Segment(polygon[i], polygon[(i+1)%nbp]);

    return oit;
  }

  template <typename OutputIterator>
  OutputIterator infinite_cvd_cell(Vertex_handle ,
                                   OutputIterator oit) const
  {
    //TODO
    return oit;
  }

  //returns true iff generators's cell is on the convex hull
  bool bvd_cell_is_infinite(const Vertex_handle generator) const
  {
    Face_circulator face = m_cdt.incident_faces(generator);
    Face_circulator begin = face;
    CGAL_For_all(face, begin){
      if(m_cdt.is_infinite(face))
        return true;
    }
    return false;
  }

};// class Constrained_voronoi_diagram


template<typename Tr>
Cvd_cell_2<Tr> dual(const Tr& tr,
                    const typename Tr::Vertex_handle& v)
{
  CGAL_triangulation_precondition( v != typename Tr::Vertex_handle());
  CGAL_triangulation_precondition( !tr.is_infinite(v));

  Constrained_voronoi_diagram_2<Tr> diagram(tr);
  return diagram.cvd_cell(v);
}

// dual(v) implementation for Delaunay_triangulation_2
//template<typename Tr, typename OutputIterator>
//OutputIterator
//dual(const Tr& tr,
//     typename Tr::Vertex_handle v,
//     OutputIterator oit)
//{
//  typedef Tr::Vertex_handle         Vertex_handle;
//  typedef Tr::Face_handle           Face_handle;
//  typedef Tr::Point                 Point;
//  typedef Tr::Segment               Segment;
//  typedef Tr::Geom_traits::Ray_2    Ray;
//  typedef Tr::Geom_traits::Vector_2 Vector_2;
//
//  CGAL_triangulation_precondition( v != Vertex_handle());
//  CGAL_triangulation_precondition( !tr.is_infinite(v));
//
//  // The Circulator moves ccw.
//  std::vector<Segment> segments;
//  std::vector<Ray> rays;
//  Tr::Face_circulator fc = tr.incident_faces(v), done(fc);
//  Point prev_cc;
//  bool first_ = true;
//  do
//  {
//    if(!tr.is_infinite(fc)) //finite edges (= segments)
//    {
//      if(first_)
//        prev_cc = tr.circumcenter(fc);
//      else
//      {
//        Point cc = tr.circumcenter(fc);
//        *oit++ = Segment(cc, prev_cc);
//        prev_cc = cc;
//      }
//      first_ = false;
//    }
//    else // infinite edges (= rays)
//    {
//      first_ = true;//for next segment
//      //find the one finite edge
//      for(int i = 0; i < 3; ++i)
//      {
//        if(!tr.is_infinite(fc,i))
//        {
//          Point m = CGAL::midpoint(fc->vertex(cw(i))->point(),
//                                   fc->vertex(ccw(i))->point());
//          Face_handle fn = fc->neighbor(i);
//
//          Point opp = fn->vertex(fn->index(fc))->point();
//          double dot_prod = (m-cc)*(m-opp);
//          if(dot_prod > 0.)      *oit++ = Ray(cc, m);
//          else if(dot_prod < 0.) *oit++ = Ray(cc, Vector(m, cc));
//          else //0. cc and m are the same point
//          {
//            Segment se = this->segment(Face_handle(fc,i));
//            Vector_2 normal = se.supporting_line().perpendicular(m).to_vector();
//            if(normal * (m - opp) > 0.)
//              *oit++ = Ray(cc, normal);
//            else
//              *oit++ = Ray(cc, -normal);
//          }
//        }
//      }
//    }
//  }
//  while(++fc != done);
//
//  return oit;
//}

} //namespace CGAL
#endif // CGAL_CONSTRAINED_VORONOI_DIAGRAM_2_H
