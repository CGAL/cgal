// Try Sylvain proposal  for a fully flexible design
// using the rebind trick


#include <CGAL/basic.h>
#include <cassert>
#include <iostream>
#include <fstream>

#define Cartesian Ca
#define Homogeneous Ho

#include "Triangulation_data_structure_2.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Triangulation_2.h>


CGAL_BEGIN_NAMESPACE

struct TDS_Bidon {
  typedef int  Vertex_handle;
  typedef int  Face_handle;
};



template <class Gt, class TDS = TDS_Bidon> 
class  Vertex_base : public Triangulation_vertex_base_2<Gt> 
{
public:
  typedef typename Gt::Point_2               Point;
  typedef Triangulation_vertex_base_2<Gt>    Base;

  typedef typename TDS::Face_handle              Face_handle;
  typedef typename TDS::Vertex_handle            Vertex_handle;

  template < class My_TDS>
  struct TDS_rebind { typedef Vertex_base<Gt, My_TDS> Rebound;};
    
  Vertex_base() : Base (){}
  Vertex_base(const Point & p, void * f = NULL) : Base(p,f) {}
      
  void set_vh(Vertex_handle vh) { _vh = vh;}
  void set_fh(Face_handle fh) {_fh  = fh ;}

  Vertex_handle get_vh() {return _vh ;}
  Face_handle   get_fh() {return _fh ;}

private:
  Vertex_handle _vh;
  Face_handle   _fh;
};


template <class Gt, class TDS = TDS_Bidon> 
class  Face_base : public Triangulation_face_base_2<Gt> 
{
public:
  typedef Triangulation_face_base_2<Gt>          Base;
 
  typedef typename TDS::Face_handle              Face_handle;
  typedef typename TDS::Vertex_handle            Vertex_handle;

  template < class My_TDS>
  struct TDS_rebind { typedef Face_base<Gt, My_TDS> Rebound; };

    
  Face_base() : Base (){}

  Face_base (void* v0, void* v1, void* v2) : Base (v0, v1, v2){}

  Face_base (void* v0, void* v1, void* v2, 
	     void* n0, void* n1, void* n2) 
    : Base (v0, v1, v2, n0, n1, n2){}
  
  void set_vh(Vertex_handle vh) { _vh = vh;}
  void set_fh(Face_handle fh) {_fh = fh;}

  Vertex_handle get_vh() {return _vh ;}
  Face_handle   get_fh() {return _fh ;}
  
private:
  Vertex_handle _vh;
  Face_handle   _fh;
};

CGAL_END_NAMESPACE

//
// A user defined Vertex 
template <class Gt, class TDS = CGAL::TDS_Bidon> 
class  My_vertex : public CGAL::Vertex_base<Gt, TDS>
{
public:
  typedef CGAL::Vertex_base<Gt,TDS>    Base;
  typedef typename TDS::Face_handle              Face_handle;
  typedef typename TDS::Vertex_handle            Vertex_handle;

  template < class My_TDS>
  struct TDS_rebind { typedef My_vertex<Gt, My_TDS> Rebound;};
    
  My_vertex() : Base (){}
  My_vertex(const Point & p, void * f = NULL) : Base(p,f) {}
      
  void set_wahou(Vertex_handle vh) { wahou = vh;}
  Vertex_handle get_wahou() {return wahou ;}

private:
 Vertex_handle wahou;
};


typedef CGAL::Cartesian<double>                         K;
typedef My_vertex<K>                                    MyVb;
typedef CGAL::Vertex_base<K>                            Vb;
typedef CGAL::Face_base<K>                              Fb;
//typedef CGAL::Triangulation_data_structure_2<Vb,Fb>     Tds;
typedef CGAL::Triangulation_data_structure_2<MyVb,Fb>     Tds;
typedef CGAL::Triangulation_2<K,Tds>                    Triangulation;

typedef K::Point_2                               Point;
typedef Triangulation::Face_handle               Face_handle;
typedef Triangulation::Vertex_handle             Vertex_handle;

typedef Triangulation::Finite_faces_iterator     Face_iterator;
typedef Triangulation::Finite_vertices_iterator  Vertex_iterator;
typedef Triangulation::Finite_edges_iterator     Edge_iterator;


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
    v_iter->set_vh(vertex);
    v_iter->set_fh(face);
    vertex = v_iter;
    face  = v_iter->face();
    v_iter->set_wahou(v_iter);
  }

  for (Face_iterator f_iter = T.faces_begin(); 
       f_iter != T.faces_end(); 
       f_iter++){
    f_iter->set_vh(vertex);
    f_iter->set_fh(face);
  }
  

  for (Vertex_iterator v_it = T.finite_vertices_begin(); 
       v_it != T.finite_vertices_end(); 
       v_it++){
    std::cerr << v_it->point() << " " << v_it->get_wahou()->point()<<"\n";
  }
  return 0;
}




