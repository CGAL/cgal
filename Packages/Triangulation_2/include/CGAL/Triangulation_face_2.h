#ifndef CGAL_TRIANGULATION_FACE_2_H
#define CGAL_TRIANGULATION_FACE_2_H

#include <CGAL/Pointer.h>
#include <CGAL/Triangulation_short_names_2.h>
#include <CGAL/Triangulation_default_data_structure_2.h>


template < class Gt, class Tds >
class CGAL_Triangulation_vertex_2;

template < class Gt, class Tds >
class CGAL_Triangulation_vertex_handle_2;

template < class Gt, class Tds >
class CGAL_Triangulation_face_handle_2;

template < class Gt, class Tds >
class CGAL_Triangulation_face_2  : public  Tds::Face
{
public:
  //  typedef Tds Tds;

  typedef Gt  Geom_traits;
  typedef typename Geom_traits::Point Point;
  typedef typename Geom_traits::Segment Segment;
  typedef typename Geom_traits::Triangle Triangle;

  typedef typename Tds::Vertex Ve;
  typedef typename Tds::Face Fa;

  typedef CGAL_Triangulation_vertex_2<Gt,Tds> Vertex;
  typedef CGAL_Triangulation_face_2<Gt,Tds> Face;

  typedef CGAL_Triangulation_vertex_handle_2<Gt,Tds> Vertex_handle;
  typedef CGAL_Triangulation_face_handle_2<Gt,Tds> Face_handle;
  //  typedef pair<Face_handle, int>     Edge;


  inline
  CGAL_Triangulation_face_2()
    : Fa()
  { }

//   inline 
//   CGAL_Triangulation_face_2(Fa f)
//     : Fa(f)
//   {}
// 
  inline
  CGAL_Triangulation_face_2(const Vertex_handle& v0,
			  const Vertex_handle& v1,
			  const Vertex_handle& v2)
    : Fa(&(*v0), &(*v1), &(*v2))
  {}
        
    
  inline
  CGAL_Triangulation_face_2(const Vertex_handle& v0,
			  const Vertex_handle& v1,
			  const Vertex_handle& v2,
			  const Face_handle& n0,
			  const Face_handle& n1,
			  const Face_handle& n2)
    : Fa(&(*v0), &(*v1), &(*v2),&(*n0), &(*n1), &(*n2)) 
  {}
        
    

  // Vertex access functions
  inline Vertex_handle vertex(int i) const
  {
    return  ((Vertex *)(Fa::vertex(i)));
  }
    
    
  inline bool has_vertex(const Vertex_handle& v) const
  {
        return (Fa::has_vertex( & (*v)) );
  }
    
    
  inline bool has_vertex(const Vertex_handle& v, int& i) const
  {
    return Fa::has_vertex( &(*v), i);
  }

  inline int index(const Vertex_handle& v) const
  {
    return Fa::index( &(*v));
  }
  
  
 

  //ACCESS FUNCTIONS
  inline
  Face_handle neighbor(int i)
  {
    return (Face *)(Fa::neighbor(i));
  }

   inline int index(Face_handle f) const
  {
    return Fa::index( &(*f));
  }
  
  inline bool has_neighbor(Face_handle f)
  {
    return Fa::has_neighbor( &(*f));
  }

  inline bool has_neighbor(Face_handle f, int i)
  {
    return Fa::has_neighbor( &(*f), i);
  }

 
  inline Face_handle handle() const
  {
    return Face_handle(this);
  }

 //Setting
  void set_vertices(Vertex_handle v0,
		    Vertex_handle v1,
		    Vertex_handle v2)
    {
        Fa::set_vertices(&(*v0), &(*v1), &(*v2));
    }
    
    void set_neighbors(Face_handle n0,
                       Face_handle n1,
                       Face_handle n2)
    {
        Fa::set_neighbor(&(*n0), &(*n1), &(*n2));
    }
    
  //void set_vertices() inherited
  //void set_neighbors() inherited

    void set_vertex(int i, Vertex_handle v)
    {
        Fa::set_vertex(i, &(*v));
    }
    
    void set_neighbor(int i, Face_handle n)
    {
        Fa::set_neighbor(i, &(*n));
    }
    
  
  //MODIFIERS
//   void insert_in_face(const Vertex_handle& v)
//   {
//     Fa::insert_in_face(&(*v));
//   }
// 
//   void insert_on_edge(const Vertex* v, int i)
//   {
//     Fa::insert_on_edge(&(*v),i);
//   }
// 
//   void insert_out(const Vertex* v, int i)
//   {
//     Fa::insert_out(&(*v),i);
//   }

//   // for compatibility with previous versions
//   void insert(const Vertex_handle& v)
//   {
//     Fa::insert_in_face(&(*v));
//   }
// 
//   // for compatibility with previous versions
//   void insert(const Vertex* v, int i)
//   {
//     Fa::insert_out(&(*v),i);
//   }
// 
//   bool remove(Vertex_handle v)
//   {
//     Fa::remove(&(*v));
//   }
// 
//   //void flip(int i) inherited
// 
//   CHECKING
//   bool is_valid(bool verbose = false, int level = 0) const inherited
// 
};

#endif CGAL_TRIANGULATION_FACE_2_H
