#ifndef CGAL_TRIANGULATION_VERTEX_2_H
#define CGAL_TRIANGULATION_VERTEX_2_H

#include <CGAL/Pointer.h>
#include <CGAL/Triangulation_default_data_structure_2.h>
#include <CGAL/Triangulation_circulators_2.h>


template < class Tds >
class CGAL_Triangulation_face_2;

template < class Tds >
class CGAL_Triangulation_vertex_handle_2;

template < class Tds >
class CGAL_Triangulation_face__handle_2;

template<class Tds>
class CGAL_Triangulation_face_circulator_2;

template<class Tds>
class CGAL_Triangulation_vertex_circulator_2;

template<class Tds>
class CGAL_Triangulation_edge_circulator_2;


template<class Tds >
class CGAL_Triangulation_vertex_2
  : public Tds::Vertex
{
public:
  
  typedef Tds Triangulation_data_structure;

  typedef typename Tds::Geom_traits  Geom_traits;
  typedef typename Geom_traits::Point Point;
  typedef typename Geom_traits::Segment Segment;
  typedef typename Geom_traits::Triangle Triangle;

  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;

  typedef CGAL_Triangulation_face_2<Tds> Face;
  typedef CGAL_Triangulation_vertex_2<Tds> Vertex;
  
  typedef CGAL_Triangulation_vertex_handle_2<Tds> Vertex_handle;
  typedef CGAL_Triangulation_face_handle_2<Tds> Face_handle;
  typedef pair<Face_handle, int>     Edge;

  typedef CGAL_Triangulation_face_circulator_2<Tds>      Face_circulator;
  typedef CGAL_Triangulation_edge_circulator_2<Tds>      Edge_circulator;
  typedef CGAL_Triangulation_vertex_circulator_2<Tds>    Vertex_circulator;

  
  CGAL_Triangulation_vertex_2()
     : Ve()
  {}

  CGAL_Triangulation_vertex_2(const Point & p)
    :  Ve(p)
  {}
    
  CGAL_Triangulation_vertex_2(const Point & p, Face_handle f)
    :  Ve(p, &(*f))
  {}

  inline void set_face(const Face_handle& f)
  {
    Ve::set_face(&(*f));
  }
    
    
  inline Face_handle face() const
  {
        return  (Face *)Ve::face();
  }
        
  inline Vertex_handle handle() const
  {
    return Vertex_handle(this);
  }
      
  //inherited :
  // cw(i), ccw(i)
  // degree(), point(), is_on_boundary, set_point(p)

  inline Vertex_circulator incident_vertices()
  {
    return Vertex_circulator(handle(), face());
  }
    
  inline 
  CGAL_Triangulation_face_circulator_2<Tds> incident_faces()
  {
    return Face_circulator(handle(), face());
  }
    
  inline Face_circulator incident_faces(const Face_handle& f)
  {
    return Face_circulator(handle(), f);
  }
    
  inline Edge_circulator incident_edges()
  {
    return Edge_circulator(handle(), face());
  }
 
  inline Edge_circulator incident_edges(const Face_handle& f)
  {
    return Edge_circulator(handle(), f);
  }
    
};
#endif CGAL_TRIANGULATION_2_H
