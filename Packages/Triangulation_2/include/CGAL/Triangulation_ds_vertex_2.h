#ifndef CGAL_TRIANGULATION_DS_VERTEX_2_H
#define CGAL_TRIANGULATION_DS_VERTEX_2_H

#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_ds_circulators_2.h>

// template <class Vb, class Fb >
// class  CGAL_Triangulation_ds_face_2 ;
// 
// template <class Vertex, class Face>
// class CGAL_Triangulation_ds_face_circulator_2;
// 
// template <class Vertex, class Face>
// class CGAL_Triangulation_ds_vertex_circulator_2;
// 
// template <class Vertex, class Face>
// class CGAL_Triangulation_ds_edge_circulator_2;
// 
template <class Vb, class Fb >
class  CGAL_Triangulation_ds_vertex_2 
  : public Vb
{
public:
  typedef typename Vb::Point Point;
  typedef CGAL_Triangulation_ds_vertex_2<Vb,Fb> Vertex;
  typedef CGAL_Triangulation_ds_face_2<Vb,Fb> Face;
  typedef pair< Face*,int> Edge;
  typedef CGAL_Triangulation_ds_face_circulator_2<Vertex,Face> Face_circulator;
  typedef CGAL_Triangulation_ds_vertex_circulator_2<Vertex,Face> Vertex_circulator;
  typedef CGAL_Triangulation_ds_edge_circulator_2<Vertex,Face> Edge_circulator;

  CGAL_Triangulation_ds_vertex_2()
    : Vb()
  {}
    
  CGAL_Triangulation_ds_vertex_2(const Point & p)
    :  Vb(p)
  {}
    
  CGAL_Triangulation_ds_vertex_2(const Point & p, Face * f)
    :  Vb(p, f )
  {}

  // set_point()
  // point()
  //inherited from Vb

  inline 
  void set_face(Face* f)
  {
    Vb::set_face(f);
  }

  inline Face* face() const
  {
    return ( (Face *) (Vb::face()) );
  }
    
  static int ccw(int i)
  {
    return (i+1) % 3;
  }
    
  static int cw(int i)
  {
    return (i+2) % 3;
  }

  // Next member function
  //assume that if the vertex is on the boundary
  //face() is the first face adjacent to the vertex
  // in counterclockwise order

//   bool is_on_boundary()
//   {
//     Face * f= face();
//     int i= f->index(this);
//     return (f->neighboor(cw(i)) == NULL);
//   }
  
  int degree() const
  {
    Face* f = face();
    
    if (f == NULL) {
      return 0;
    }
    int i = f->index(this);
    
    Face* ptr1 = f->neighbor(ccw(i));
    Face* ptr2 = f;
    f = f->neighbor(cw(i));
    
    int count = 2;
    while(ptr1 != f){
      count++;
      i = ptr1->index(ptr2);
      ptr2 = ptr1;
      ptr1 = ptr1->neighbor(cw(i));
    }
    return count;
  }
  
  inline Vertex_circulator incident_vertices()
  {
    return Vertex_circulator(this, face());
  }
    
  inline Face_circulator incident_faces()
  {
    return Face_circulator(this, face());
  }
    
  inline Face_circulator incident_faces(const Face* f)
  {
    return Face_circulator(this, f);
  }
    
  inline Edge_circulator incident_edges()
  {
    return Edge_circulator(this, face());
  }
       
  inline Edge_circulator incident_edges(const Face* f)
  {
    return Edge_circulator(this, f);
  }
    
   bool is_valid(bool verbose = false, int level = 0) const
  {
    bool result = Vb::is_valid();
    result = result && face()->has_vertex(this);
    return result;
  }
};

#endif CGAL_TRIANGULATION_DS_VERTEX_2_H
