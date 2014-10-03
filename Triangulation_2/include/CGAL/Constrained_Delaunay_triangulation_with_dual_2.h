#ifndef CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_DUAL_2_H
#define CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_DUAL_2_H

#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Delaunay_triangulation_vertex_base_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_face_base_2.h>

#include <CGAL/Polygon_2.h>

#include <utility>

namespace CGAL {

template <class Cdt>
class Bvd_cell : 
  public std::pair<typename Cdt::Vertex_handle, // generator
                   typename Cdt::Polygon> // cell
{
  typedef std::pair<typename Cdt::Vertex_handle,
                    typename Cdt::Polygon>           Base;
  
  typedef typename Cdt::Vertex_handle    Vertex_handle;
public:
  typedef typename Cdt::Polygon          Polygon;

public:
  Bvd_cell()
    : Base()
  {
  }
  Bvd_cell(Vertex_handle v, const Polygon& poly)
    : Base(v, poly)
  {
  }
                  
  bool operator<(const Bvd_cell& cell) const
  {
    return m_vertex->point() < cell.vertex()->point();
  };

  Vertex_handle vertex() const
  {
    return this->first;
  }

  typename Cdt::Polygon cell() const
  {
    return this->second;
  }
}; //end CLASS Bvd_cell


template <class Gt,
          class Tds = Triangulation_data_structure_2<
                        Delaunay_triangulation_vertex_base_2<Gt>,
                        Constrained_Delaunay_triangulation_face_base_2<Gt> >,
      class Itag = No_intersection_tag >
class Constrained_Delaunay_triangulation_with_dual_2
  : public Constrained_Delaunay_triangulation_2<Gt, Tds, Itag>
    //public CGAL::Constrained_triangulation_plus_2<Cdt> 
{
public:
  // typedefs for basic primitives 
  typedef Constrained_Delaunay_triangulation_2<Gt,Tds,Itag>    Base;
  typedef Constrained_Delaunay_triangulation_with_dual_2 Cdt;
  typedef typename Base::Geom_traits       Geom_traits;
  typedef typename Base::Intersection_tag Intersection_tag;

  typedef typename Base::Constraint             Constraint;
  typedef typename Base::Vertex_handle          Vertex_handle;
  typedef typename Base::Face_handle            Face_handle;
  typedef typename Base::Edge                   Edge;
  typedef typename Base::Finite_faces_iterator  Finite_faces_iterator;
  typedef typename Base::Face_circulator        Face_circulator;
  typedef typename Base::Edge_circulator        Edge_circulator;
  typedef typename Base::size_type              size_type;
  typedef typename Base::Locate_type            Locate_type;

  typedef typename Cdt::Geom_traits             Kernel;
  typedef typename Kernel::FT                   FT;
  typedef typename Kernel::Point_2              Point;
  typedef typename Kernel::Vector_2             Vector;
  typedef typename Kernel::Line_2               Line;
  typedef typename Cdt::Segment                 Segment;
  typedef typename Cdt::Triangle                Triangle;

  // typedefs for Bounded Voronoi Diagram (Bvd)
  typedef Bvd_cell<Cdt> Bvd_cell;
  typedef typename Bvd_cell::Polygon Polygon;
  typedef typename std::list<Bvd_cell> Bvd; 

protected:
    FT m_bounding_box[4]; // xmin,xmax,ymin,ymax
    Bvd m_bvd;

public:
    Cdt() 
    {
      m_bounding_box[0] = m_bounding_box[2] = 0.0;
      m_bounding_box[1] = m_bounding_box[3] = 1.0;
    }

  enum {INSIDE = -1,
        UNDETERMINED = 0,
        OUTSIDE = 1};

  //----------------------------------------------------------------
  //--------------------ABOUT FACES SIGHT---------------------------
  //----------------------------------------------------------------

public:

  // blind = false IFF each face sees its circumcenter
  void tag_all_faces_blind(const bool blind) 
  {
    for(All_faces_iterator f = this->all_faces_begin();
         f != this->all_faces_end();
         ++f)
      f->blind() = blind;
  }

  // blind test for each face
  // if true, set corresponding barrier constraint
  void tag_faces_blind()
  {
    if(dimension() < 2)
      return;

    tag_all_faces_blind(false);

    // for each constrained edge, mark blinded triangles
    for(Finite_edges_iterator e = this->finite_edges_begin();
         e != this->finite_edges_end();
         ++e)
    {
      Edge edge = *e;
      if(this->is_constrained(edge))
      {
        tag_neighbors_blind(edge);
        tag_neighbors_blind(this->mirror_edge(edge));
      }
    }
  }

  // returns the edge seen from the "other side"
  Edge mirror_edge(const Edge& e) const
  {
    Face_handle neighb = e.first->neighbor(e.second);
    return Edge(neighb, neighb->index(e.first));
  }

  // test face for blindness with respect to the edge constraint
  void tag_face_blind(Face_handle& f, const Edge& constraint)
  {  
    if(segment_hides_circumcenter(this->segment(constraint),
                                               this->triangle(f)))
    {
      f->blind() = true;
      f->blinding_constraint() = constraint;
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
    CGAL_assertion(this->is_constrained(constraint));
    Face_handle seed = constraint.first;

    if(!this->is_infinite(seed) 
       && !seed->blind() 
       && triangle(seed).area() != 0) //to avoid flat triangles outside the domain
    {
      std::stack<Face_handle> faces;
      faces.push(seed);

      while(!faces.empty())
      {
        Face_handle f = faces.top();
        faces.pop();
        this->tag_face_blind(f, constraint);
        if( f->blind())
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
      if(!this->is_constrained(edge_i) &&
          !fi->blind() &&
          !is_infinite(fi)) 
        faces.push(fi);
    }
  }

  void tag_faces_location(const int location)
  {
    for(All_faces_iterator f = all_faces_begin();
      f != all_faces_end();
      f++)
      f->location() = location;
  }


  /*--------------------------------------------------------------
  ---------------------- BVD CONSTRUCTION ------------------------
  --------------------------------------------------------------*/

  void construct_bvd()
  {
    m_bvd.clear();

    tag_faces_blind();

    for(Finite_vertices_iterator v = this->finite_vertices_begin();
        v != this->finite_vertices_end();
        ++v)
    {
      if(!this->bvd_cell_is_infinite(v))
        m_bvd.push_back(this->bvd_cell(v));
    }
  }


  // assemble a cell of the bounded Voronoi diagram
  // incident to vertex v
  Bvd_cell bvd_cell(Vertex_handle v) const
  {
    Polygon polygon;

    CGAL_assertion(!is_infinite(v));
    Face_circulator face = this->incident_faces(v);
    Face_circulator end = face;
    Face_circulator next = face;

    CGAL_For_all(face, end)
    {
      next++;
      Line line(this->circumcenter(face), this->circumcenter(next));
      Point intersect;

      if(!face->blind()) //face sees
      {
        polygon.push_back(this->circumcenter(face));
        if(next->blind())  //next doesn't
        {
          CGAL_assertion(do_intersect(line, this->segment(next->blinding_constraint())));
          CGAL::assign(intersect,
            CGAL::intersection(line, Line(this->segment(next->blinding_constraint()))));
          polygon.push_back(intersect);
        }
      }
      else //face doesn't see
      {
        if(!next->blind()) //next sees
        {
          CGAL_assertion(do_intersect(line, this->segment(face->blinding_constraint())));
          CGAL::assign(intersect,
            CGAL::intersection(line, Line(this->segment(face->blinding_constraint()))));
          polygon.push_back(intersect);
        }
        else //next doesn't
        {
          if(face->blinding_constraint() != next->blinding_constraint()
            && face->blinding_constraint() != this->mirror_edge(next->blinding_constraint()))
            // the 2 blinding_constraints are different
          {
            CGAL_assertion(do_intersect(line, this->segment(face->blinding_constraint())));
            CGAL::assign(intersect,
              CGAL::intersection(line, Line(this->segment(face->blinding_constraint()))));
            polygon.push_back(intersect);

            Point intersection2;
            CGAL_assertion(do_intersect(line, this->segment(next->blinding_constraint())));
            CGAL::assign(intersection2,
              CGAL::intersection(line, Line(this->segment(next->blinding_constraint()))));
            polygon.push_back(intersection2);
          }
          //else: it's the same constraint--> do nothing
        }
      }
    }//end CGAL_For_all

    return Bvd_cell(v, polygon); 
  }  


  ////returns true IFF generator's cell is clipped by at least one constrained Edge
  //bool bvd_cell_is_clipped(const Vertex_handle generator)
  //{
  //  Face_circulator face = this->incident_faces(generator);
  //  Face_circulator begin = face;
  //  CGAL_For_all(face, begin){
  //    if(face->blind())
  //      return true;
  //  }
  //  return false;
  //}

  //returns true IFF generators's cell is finite
  bool bvd_cell_is_infinite(const Vertex_handle generator)
  {
    Face_circulator face = this->incident_faces(generator);
    Face_circulator begin = face;
    CGAL_For_all(face, begin){
      if(this->is_infinite(face))
        return true;
    }
    return false;
  }


//----------------------------------------------------
//------------------ TOOLS ---------------------------
//----------------------------------------------------

  bool is_inside(const Point& query)
  {
    Face_handle f = locate(query);
    if(f == NULL)
      return false;
    if(f->location() == INSIDE)
      return true;
    return false;
  }

  bool incident_constraints(Vertex_handle v)
  {
    Edge_circulator e = incident_edges(v);
    Edge_circulator end = e;
    CGAL_For_all(e,end)
      if((*e).first->is_constrained((*e).second))
        return true;
    return false;
  }

  //void tag_faces_inside_outside(std::list<Point>& seeds)
  //{
  //  // reset all tags undetermined
  //  tag_faces_location(INSIDE);

  //  // pick infinite seed face
  //  Face_handle seed = infinite_vertex()->face();
  //  seed->location() = OUTSIDE;
  //  std::stack<Face_handle> faces;
  //  faces.push(seed);

  //  // locate from all seeds
  //  std::list<Point>::iterator it;
  //  for(it = seeds.begin();
  //    it != seeds.end();
  //    it++)
  //  {
  //    const Point& p = *it;
  //    Face_handle f = locate(p);
  //    if(f != NULL)
  //      faces.push(f);
  //  }

  //  while(!faces.empty())
  //  {
  //    Face_handle f = faces.top();
  //    faces.pop();
  //    const int& location = f->location();
  //    for(unsigned int i=0;i<3;i++)
  //    {
  //      if(f->neighbor(i) != NULL)
  //        if(f->neighbor(i)->location() == INSIDE && 
  //          !f->is_constrained(i))
  //        {
  //          f->neighbor(i)->location() = OUTSIDE;
  //          faces.push(f->neighbor(i));
  //        }
  //    }
  //  }
  //}

  //void tag_faces_inside_outside()
  //{
  //  // reset all tags undetermined
  //  tag_faces_location(UNDETERMINED);

  //  // pick one seed face
  //  Face_handle seed = infinite_vertex()->face();
  //  seed->location() = OUTSIDE;
  //  std::stack<Face_handle> faces;
  //  faces.push(seed);

  //  while(!faces.empty())
  //  {
  //    Face_handle f = faces.top();
  //    faces.pop();
  //    const int& location = f->location();
  //    ASSERT(location == OUTSIDE || location == INSIDE);
  //    for(unsigned int i=0;i<3;i++)
  //    {
  //      if(f->neighbor(i) != NULL)
  //        if(f->neighbor(i)->location() == UNDETERMINED)
  //        {
  //          // swap inside / outside if crosses a tagged facet
  //          bool constrained = f->is_constrained(i);
  //          f->neighbor(i)->location() = constrained ? -location : location;
  //          faces.push(f->neighbor(i));
  //        }
  //    }
  //  }
  //}

  // compute bounding box
  void compute_bounding_box()
   {
      Finite_vertices_iterator v;
      for(v = this->finite_vertices_begin();
          v != this->finite_vertices_end();
          v++)
      {
        if(v == this->finite_vertices_begin())
        {
          m_bounding_box[0] = m_bounding_box[1] = v->point().x();
          m_bounding_box[2] = m_bounding_box[3] = v->point().y();
        }
        else
        {
          m_bounding_box[0] = std::min(m_bounding_box[0],v->point().x());
          m_bounding_box[1] = std::max(m_bounding_box[1],v->point().x());
          m_bounding_box[2] = std::min(m_bounding_box[2],v->point().y());
          m_bounding_box[3] = std::max(m_bounding_box[3],v->point().y());
        }
      }
      CGAL_assertion(valid_bounding_box());
   }

   bool valid_bounding_box() const
   {
     return (m_bounding_box[1] - m_bounding_box[0]) > 0.0 &&
            (m_bounding_box[3] - m_bounding_box[2]) > 0.0;
   }

};// class CCDT

} //namespace CGAL
#endif // CGAL_CONSTRAINED_DELAUNAY_TRIANGULATION_WITH_DUAL_2_H