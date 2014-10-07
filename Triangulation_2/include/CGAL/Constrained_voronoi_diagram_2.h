#ifndef CGAL_CONSTRAINED_VORONOI_DIAGRAM_2_H
#define CGAL_CONSTRAINED_VORONOI_DIAGRAM_2_H

#include <CGAL/Polygon_2.h>

#include <utility>

namespace CGAL {

template <class Cdt>
class Cvd_cell
{
  typedef typename Cdt::Vertex_handle    Vertex_handle;

public:
  typedef typename Cdt::Geom_traits::Segment_2    Segment;
  typedef typename Cdt::Geom_traits::Ray_2        Ray;
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
  Cvd_cell(Vertex_handle v)
    : m_vertex(v)
    , m_segments()
    , m_rays(2)
    , m_is_valid(false)
  { }

  bool operator<(const Cvd_cell& cell) const
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
  //access iterators
  typedef typename std::vector<Segment>::iterator       segment_iterator;
  typedef typename std::vector<Ray>::iterator           ray_iterator;

  segment_iterator segments_begin()       { return m_segments.begin(); }
  segment_iterator segments_end()         { return m_segments.end(); }

  ray_iterator rays_begin()       { return m_rays.begin(); }
  ray_iterator rays_end()         { return m_rays.end(); }

}; //end CLASS Cvd_cell


// Cdt should be of the type Constrained_Delaunay_triangulation_face_base_2
template <class Cdt>
class Constrained_voronoi_diagram_2
  : public std::list< Cvd_cell<Cdt> >
{
  typedef std::list< Cvd_cell<Cdt> > Base;

public:
  typedef Constrained_voronoi_diagram_2<Cdt>         Cvd; 
  typedef typename Cvd_cell<Cdt>                     Cvd_cell;
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

  typedef typename Cdt::Geom_traits             Kernel;
  typedef typename Kernel::FT                   FT;
  typedef typename Kernel::Point_2              Point;
  typedef typename Kernel::Vector_2             Vector;
  typedef typename Kernel::Line_2               Line;
  typedef typename Cdt::Segment                 Segment;
  typedef typename Cdt::Triangle                Triangle;

  typedef typename Base::iterator       iterator;
  typedef typename Base::const_iterator const_iterator;

protected:
  const Cdt* m_pCdt;

public:
  Cvd(const Cdt* p_cdt)
    : Base()
    , m_pCdt(p_cdt)
  {
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
    for(All_faces_iterator f = m_pCdt->all_faces_begin();
         f != m_pCdt->all_faces_end();
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
    for(Finite_edges_iterator e = m_pCdt->finite_edges_begin();
         e != m_pCdt->finite_edges_end();
         ++e)
    {
      Edge edge = *e;
      if(m_pCdt->is_constrained(edge))
      {
        tag_neighbors_blind(edge);
        tag_neighbors_blind(m_pCdt->mirror_edge(edge));
      }
    }
  }

  // test face for blindness with respect to the edge constraint
  void tag_face_blind(Face_handle& f, const Edge& constraint)
  {  
    if(segment_hides_circumcenter(m_pCdt->segment(constraint),
                                  m_pCdt->triangle(f)))
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
    CGAL_assertion(m_pCdt->is_constrained(constraint));
    Face_handle seed = constraint.first;

    if(!m_pCdt->is_infinite(seed) 
       && !seed->blind() 
       && triangle(seed).area() != 0) //to avoid flat triangles outside the domain
    {
      std::stack<Face_handle> faces;
      faces.push(seed);

      while(!faces.empty())
      {
        Face_handle f = faces.top();
        faces.pop();
        m_pCdt->tag_face_blind(f, constraint);
        if( f->blind())
          m_pCdt->push_unvisited_neighbors(f, faces);
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
      if(!m_pCdt->is_constrained(edge_i) &&
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
    this->clear();
    tag_faces_blind();
        
    for(Finite_vertices_iterator v = m_pCdt->finite_vertices_begin();
        v != m_pCdt->finite_vertices_end();
        ++v)
    {
      this->push_back(m_pCdt->cvd_cell(v));
    }
  }

  // assemble a cell of the bounded Voronoi diagram
  // incident to vertex v
public:
  template<typename OutputIterator>
  OutputIterator cvd_cell(Vertex_handle v, OutputIterator oit) const
  {
    if(bvd_cell_is_infinite(v))
      infinite_cvd_cell(v, oit);
    else
      finite_cvd_cell(v, oit);
    return oit;
  }

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

private:
  template <typename OutputIterator>
  OutputIterator finite_cvd_cell(Vertex_handle v, OutputIterator oit) const
  {
    std::vector<Point> polygon;
    
    CGAL_assertion(!m_pCdt->is_infinite(v));
    Face_circulator face = m_pCdt->incident_faces(v);
    Face_circulator end = face;
    Face_circulator next = face;

    CGAL_For_all(face, end)
    {
      next++;
      Line line(m_pCdt->circumcenter(face), m_pCdt->circumcenter(next));
      Point intersect;

      if(!face->blind()) //face sees
      {
        polygon.push_back(m_pCdt->circumcenter(face));
        if(next->blind())  //next doesn't
        {
          CGAL_assertion(do_intersect(line, m_pCdt->segment(next->blinding_constraint())));
          CGAL::assign(intersect,
            CGAL::intersection(line, Line(m_pCdt->segment(next->blinding_constraint()))));
          polygon.push_back(intersect);
        }
      }
      else //face doesn't see
      {
        if(!next->blind()) //next sees
        {
          CGAL_assertion(do_intersect(line, m_pCdt->segment(face->blinding_constraint())));
          CGAL::assign(intersect,
            CGAL::intersection(line, Line(m_pCdt->segment(face->blinding_constraint()))));
          polygon.push_back(intersect);
        }
        else //next doesn't
        {
          if(face->blinding_constraint() != next->blinding_constraint()
            && face->blinding_constraint() != m_pCdt->mirror_edge(next->blinding_constraint()))
            // the 2 blinding_constraints are different
          {
            CGAL_assertion(do_intersect(line, m_pCdt->segment(face->blinding_constraint())));
            CGAL::assign(intersect,
              CGAL::intersection(line, Line(m_pCdt->segment(face->blinding_constraint()))));
            polygon.push_back(intersect);

            Point intersection2;
            CGAL_assertion(do_intersect(line, m_pCdt->segment(next->blinding_constraint())));
            CGAL::assign(intersection2,
              CGAL::intersection(line, Line(m_pCdt->segment(next->blinding_constraint()))));
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
  OutputIterator infinite_cvd_cell(Vertex_handle v,
                                   OutputIterator oit) const
  {
    return oit;
  }

  ////returns true IFF generator's cell is clipped by at least one constrained Edge
  //bool bvd_cell_is_clipped(const Vertex_handle generator)
  //{
  //  Face_circulator face = m_pCdt->incident_faces(generator);
  //  Face_circulator begin = face;
  //  CGAL_For_all(face, begin){
  //    if(face->blind())
  //      return true;
  //  }
  //  return false;
  //}

  //returns true iff generators's cell is on the convex hull
  bool bvd_cell_is_infinite(const Vertex_handle generator) const
  {
    Face_circulator face = m_pCdt->incident_faces(generator);
    Face_circulator begin = face;
    CGAL_For_all(face, begin){
      if(m_pCdt->is_infinite(face))
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


};// class Constrained_voronoi_diagram

} //namespace CGAL
#endif // CGAL_CONSTRAINED_VORONOI_DIAGRAM_2_H