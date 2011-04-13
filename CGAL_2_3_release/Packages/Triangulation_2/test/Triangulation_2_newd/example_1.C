// Try to make deferred instantiation at the level of the Tds


#include <CGAL/basic.h>
#include <cassert>
#include <iostream>
#include <fstream>

#define Cartesian Ca
#define Homogeneous Ho

#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>

CGAL_BEGIN_NAMESPACE

template <class Gt, class Ref> 
class  My_vertex : public Triangulation_vertex_base_2<Gt> 
{
public:
  typedef typename Gt::Point_2               Point;
  typedef Triangulation_vertex_base_2<Gt>    Vertex_base;

  typedef typename Ref::Face                     Tds_Face;
  typedef typename Ref::Vertex                   Tds_Vertex;
  typedef Triangulation_face_2<Gt,Ref>           Face;
  typedef Triangulation_vertex_2<Gt,Ref>         Vertex;
  typedef Triangulation_face_handle_2<Gt,Ref>    Face_handle;
  typedef Triangulation_vertex_handle_2<Gt,Ref>  Vertex_handle;
    
  My_vertex() : Vertex_base (){}

  My_vertex (const Point & p, void * f = NULL) : Vertex_base (p, f) {}
    
  void set_v (Vertex* v) {_v = v;}
  void set_f (Face* f) {_f = f;}
  void set_vh(Vertex_handle vh) { _v = &(*vh);}
  void set_fh(Face_handle fh) {_f = &(*fh);}

  Vertex* get_v() { return _v;}
  Face*   get_f() { return _f;}
  Vertex_handle get_vh() {return Vertex_handle(_v) ;}
  Face_handle   get_fh() {return Face_handle(_f);}

private:
  Vertex*  _v;
  Face*    _f;
};


template <class Gt, class Ref> 
class  My_face : public Triangulation_face_base_2<Gt> 
{
public:
  typedef Triangulation_face_base_2<Gt>            Face_base;
 
  typedef typename Ref::Face                     Tds_Face;
  typedef typename Ref::Vertex                   Tds_Vertex;
  typedef Triangulation_face_2<Gt,Ref>           Face;
  typedef Triangulation_vertex_2<Gt,Ref>         Vertex;
  typedef Triangulation_face_handle_2<Gt,Ref>    Face_handle;
  typedef Triangulation_vertex_handle_2<Gt,Ref>  Vertex_handle;
  
  My_face() : Face_base (){}

  My_face (void* v0, void* v1, void* v2) : Face_base (v0, v1, v2){}

  My_face (void* v0, void* v1, void* v2, 
	   void* n0, void* n1, void* n2) 
    : Face_base (v0, v1, v2, n0, n1, n2){}
  
  void set_v (Vertex* v) {_v = v;}
  void set_f (Face* f) {_f = f;}
  void set_vh(Vertex_handle vh) { _v = &(*vh);}
  void set_fh(Face_handle fh) {_f = &(*fh);}

  Vertex* get_v() { return _v;}
  Face*   get_f() { return _f;}
  Vertex_handle get_vh() {return Vertex_handle(_v) ;}
  Face_handle   get_fh() {return Face_handle(_f);}
  
private:
  Vertex*  _v;
  Face*    _f;
};

//
//The user who wants to use its own face or own vertex
// has to rewrite this class
template <class Gt>
class My_TDS:
 public  Triangulation_default_data_structure_2 <
                          Gt , 
                          My_vertex <Gt, My_TDS<Gt> > ,
                          My_face <Gt, My_TDS<Gt> > 
                          >

{
public:  // CREATION
  My_TDS() {} 
};

CGAL_END_NAMESPACE

typedef double coord_type;
typedef CGAL::Cartesian<coord_type>  Rpst;
typedef CGAL::Triangulation_euclidean_traits_2<Rpst> Gt;
typedef CGAL::My_TDS<Gt>                             My_tds;
typedef CGAL::Triangulation_2<Gt,My_tds>             Triangulation;

typedef Gt::Point_2          Point;
typedef Triangulation::Face  Face;
typedef Triangulation::Vertex Vertex;
typedef Triangulation::Face_handle  Face_handle;
typedef Triangulation::Vertex_handle Vertex_handle;

typedef Triangulation::Face_circulator  Face_circulator;
typedef Triangulation::Vertex_circulator  Vertex_circulator;

typedef Triangulation::Locate_type Locate_type;

typedef Triangulation::Face_iterator  Face_iterator;
typedef Triangulation::Vertex_iterator  Vertex_iterator;
typedef Triangulation::Edge_iterator  Edge_iterator;
typedef Triangulation::Line_face_circulator  Line_face_circulator;


int  main()
{
  Triangulation T;
  std::vector<Point> points;

  points.push_back(Point(0,0));
  points.push_back(Point(1,1));
  points.push_back(Point(2,2));

  for (unsigned int i = 0; i < points.size(); i++)
    T.insert(points[i]);
  
  Vertex_iterator v0_iter = T.finite_vertices_begin();
  
  Vertex_handle vertex = T.finite_vertex();
  Face_handle face = v0_iter->face();
  
  for (Vertex_iterator v_iter = T.finite_vertices_begin(); 
       v_iter != T.finite_vertices_end(); 
       v_iter++){
    v_iter->set_v(&(*vertex));
    v_iter->set_vh(vertex);
    v_iter->set_f(&(*face));
    v_iter->set_fh(face);
  }

  for (Face_iterator f_iter = T.faces_begin(); 
       f_iter != T.faces_end(); 
       f_iter++){
    f_iter->set_v(&(*vertex));
    f_iter->set_vh(vertex);
    f_iter->set_f(&(*face));
    f_iter->set_fh(face);
  }
  

  Vertex_iterator v_it = T.finite_vertices_begin();
  std::cout<< v_it->get_v()->point()<<"\n";
  std::cout<< v_it->get_vh()->point()<<"\n";

  return 0;
}




