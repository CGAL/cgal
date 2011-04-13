#include <fstream>
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Point_2.h>
#ifndef _MSC_VER 
#include <CGAL/Nef_2/HDS_items.h>
#include <CGAL/Nef_2/HalfedgeDS_using_in_place_list.h>
#else
#include <CGAL/Nef_2/MSVC_HalfedgeDS_using_in_place_list.h>
#endif
#include <CGAL/Nef_2/PM_decorator.h>

struct HDS_traits { 
  typedef CGAL::Cartesian<double> Kernel;
  typedef Kernel::Point_2 Point; 
  typedef bool Mark;
};

typedef  HalfedgeDS_using_in_place_list<HDS_traits,HDS_items> HDS;
typedef  CGAL::PM_decorator< HDS > PM_dec;
typedef  PM_dec::Point           Point;
typedef  PM_dec::Vertex_handle   v_handle;
typedef  PM_dec::Halfedge_handle e_handle;
typedef  PM_dec::Face_handle     f_handle;
typedef  PM_dec::Vertex_base     Vertex_base;
typedef  PM_dec::Point_const_iterator Point_const_iterator;

int main()
{
  SETDTHREAD(13);
  HDS  H;
  PM_dec D(H);
  v_handle v1 = D.new_vertex(Vertex_base(Point(30,30)));
  v_handle v2 = D.new_vertex(Vertex_base(Point(30,230)));
  v_handle v3 = D.new_vertex(Vertex_base(Point(230,230)));
  v_handle v4 = D.new_vertex(Vertex_base(Point(230,30)));

  e_handle e1 = D.new_halfedge_pair(v1,v2);
  e_handle e2 = D.new_halfedge_pair(v3,v4);
                D.new_halfedge_pair(e2,D.twin(e1));
                D.new_halfedge_pair(v1,v4);

  Point_const_iterator pit = D.points_begin();
  Point p = *pit;
  D.print_statistics();
  D.check_integrity_and_topological_planarity(false);
  D.delete_vertex(v2);
  D.check_integrity_and_topological_planarity(false);
  D.print_statistics();
  CGAL::print_as_leda_graph(std::cerr,D,CGAL::KERNELPNT());
  return 0;
}


