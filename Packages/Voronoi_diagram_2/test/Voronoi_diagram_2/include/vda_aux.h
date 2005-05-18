#ifndef VDA_AUX_H
#define VDA_AUX_H 1

#include <CGAL/basic.h>
#include <iostream>
#include <CGAL/Triangulation_utils_2.h>

template<class Vertex_handle, class Site_t>
struct Project_site
{
  typedef Site_t Site_2;
  Site_2 operator()(const Vertex_handle& v) const
  {
    return v->site();
  }
};


template<class Vertex_handle, class Site_t>
struct Project_point
{
  typedef Site_t Site_2;
  Site_2 operator()(const Vertex_handle& v) const
  {
    return v->point();
  }
};


template<class VDA, class Point_t>
struct Project_dual
{
  typedef Point_t Point_2;
  typedef typename VDA::Dual_face_handle Face_handle;

  Point_2 operator()(const VDA& vda, const Face_handle& f) const
  {
    return vda.dual().dual( f );
  }
};

template<class VDA, class Point_t>
struct Project_primal
{
  typedef Point_t Point_2;
  typedef typename VDA::Dual_face_handle Face_handle;

  Point_2 operator()(const VDA& vda, const Face_handle& f) const
  {
    return vda.dual().primal( f );
  }
};

template<class VDA, class Site_t>
struct Project_ag_dual
{
  typedef Site_t Site_2;
  typedef typename Site_2::Point_2 Point_2;
  typedef typename VDA::Dual_face_handle Face_handle;

  Point_2 operator()(const VDA& vda, const Face_handle& f) const
  {
    CGAL::Object o = vda.dual().dual(f);
    Site_2 s;
    if ( CGAL::assign(s, o) ) {
      return s.point();
    } else{
      assert( false );
      return Point_2();
    }
  }
};





template<class VDA, class Projector>
void print_halfedge(const VDA& vda,
		    const typename VDA::Halfedge_handle& he,
		    const Projector& project,
		    std::ostream& os = std::cout) {
  typename VDA::Dual_edge e = he->dual_edge();
  print_dual_edge(vda, e, project, os);
}

template<class VDA, class Projector>
void print_halfedge(const VDA& vda,
		    const typename VDA::Halfedge& he,
		    const Projector& project,
		    std::ostream& os = std::cout) {
  typename VDA::Dual_edge e = he.dual_edge();
  print_dual_edge(vda, e, project, os);
}


template<class VDA, class Projector>
void print_dual_edge(const VDA& vda,
		     const typename VDA::Dual_edge& e,
		     const Projector& project,
		     std::ostream& os = std::cout) {
  typedef CGAL::Triangulation_cw_ccw_2  CW_CCW_2;

  if ( vda.dual().is_infinite( e.first->vertex( CW_CCW_2::cw(e.second) ) )
       ) {
    os << "inf - " << std::flush;
  } else {
    os << project(  e.first->vertex( CW_CCW_2::cw(e.second) )  )
       << " - " << std::flush;
  }
  if ( vda.dual().is_infinite( e.first->vertex( CW_CCW_2::ccw(e.second) ) )
       ) {
    os << "inf" << std::endl;
  } else {
    os << project(  e.first->vertex( CW_CCW_2::ccw(e.second) )  )
       << std::endl;
  }
}

template<class VDA, class Projector>
void print_dual_vertex(const VDA& vda,
		       const typename VDA::Dual_vertex_handle& v,
		       const Projector& project,
		       std::ostream& os = std::cout) {
   if ( vda.dual().is_infinite(v) ) {
     os << "inf" << std::endl;
   } else {
     os << project(v) << std::endl;
   }
}

#endif
