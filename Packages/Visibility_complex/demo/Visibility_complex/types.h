#ifndef VISIBILITY_COMPLEX_ARR_TYPES_H
#define VISIBILITY_COMPLEX_ARR_TYPES_H

// kernel definition
typedef double FT;
typedef CGAL::Simple_cartesian<FT> K1;
typedef CGAL::Filtered_kernel<K1> Kernel;
struct K : public Kernel {};

// needed kernel objects
typedef K::Point_2 Point_2;
typedef K::Segment_2 Segment_2;
typedef K::Line_2 Line_2;

// visibility complex definition
typedef CGAL::Visibility_complex_point_traits<K> Gt;
typedef CGAL::Visibility_complex_2<Gt> Visibility_complex;
typedef Visibility_complex::Antichain Antichain;
typedef Antichain::Cw_traits Cw_traits;
typedef Visibility_complex::Vertex_handle Vertex_handle;
typedef Visibility_complex::Vertex Vertex;
typedef Visibility_complex::Face_handle Face_handle;
typedef Visibility_complex::Face Face;
typedef Visibility_complex::Edge_handle Edge_handle;
typedef Visibility_complex::Edge Edge;

// special vertex
template < class Gt, class Vb = CGAL::Triangulation_vertex_base_2<Gt> >
class My_vertex_base 
  : public  Vb
{
  typedef Vb                              Base;
public:
  typedef typename Vb::Face_handle        DT_Face_handle;
  typedef typename Vb::Point              Point;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS2>::Other    Vb2;
    typedef My_vertex_base<Gt,Vb2>                           Other;
  };

private:
  Vertex_handle  va_;

public:
  My_vertex_base() : Base() {}
  My_vertex_base(const Point & p) : Base(p) {}
  My_vertex_base(const Point & p, DT_Face_handle f) : Base(f,p) {}
  My_vertex_base(DT_Face_handle f) : Base(f) {}

  void set_associated_vertex(Vertex_handle va) { va_ = va;}
  Vertex_handle get_associated_vertex() {return va_ ; }
};

// Delaunay triangulations definition
typedef CGAL::Delaunay_triangulation_2<K> DT;

typedef CGAL::Triangulation_data_structure_2<My_vertex_base<K> > Special_tds;
typedef CGAL::Delaunay_triangulation_2<K, Special_tds> Special_DT;
typedef DT::Vertex_handle DT_vertex_handle;
typedef Special_DT::Vertex_handle Special_DT_vertex_handle;
typedef DT::Vertex DT_vertex;
typedef Special_DT::Vertex Special_DT_vertex;

#endif
