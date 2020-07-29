#define CGAL_CDT_2_DEBUG_INTERSECTIONS 1
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>

typedef CGAL::Epick Kernel;
typedef Kernel::FT FieldNumberType;
typedef Kernel::Point_2 Point2;
typedef Kernel::Point_3 Point3;

template <typename Vb>
class My_vertex_base : public Vb {
  std::size_t time_stamp_;
public:
  My_vertex_base() : Vb(), time_stamp_(-1) {
  }

  My_vertex_base(const My_vertex_base& other) :
    Vb(other),
    time_stamp_(other.time_stamp_)
  {}

  typedef CGAL::Tag_true Has_timestamp;

  std::size_t time_stamp() const {
    return time_stamp_;
  }
  void set_time_stamp(const std::size_t& ts) {
    time_stamp_ = ts;
  }

  template < class TDS >
  struct Rebind_TDS {
    typedef typename Vb::template Rebind_TDS<TDS>::Other Vb2;
    typedef My_vertex_base<Vb2> Other;
  };
};

struct FaceInfo2
{
  unsigned long long m_id;
};
typedef CGAL::Projection_traits_xy_3<Kernel> TriangulationTraits;
typedef CGAL::Triangulation_vertex_base_with_id_2<TriangulationTraits> VertexBaseWithId;
typedef My_vertex_base<VertexBaseWithId> Vb2;
typedef CGAL::Triangulation_vertex_base_2<TriangulationTraits, Vb2> VertexBase;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, Kernel> FaceBaseWithInfo;
typedef CGAL::Constrained_triangulation_face_base_2<TriangulationTraits, FaceBaseWithInfo> FaceBase;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBase> TriangulationData;
typedef CGAL::Constrained_Delaunay_triangulation_2<TriangulationTraits, TriangulationData, CGAL::Exact_predicates_tag> ConstrainedTriangulation;
    typedef CGAL::Constrained_triangulation_plus_2<ConstrainedTriangulation> CDT;

int main()
{
  CDT cdt;
  const Point3 A(539.5294108288881, 332.45151278002265, 109.660400390625);
  const Point3 B(538.779296875, 329.10546875, 109.707275390625);
  const Point3 C(539.74609375, 332.431640625, 109.660400390625);
  const Point3 D(539.68266371486789, 333.13513011783971, 109.649658203125);
  const Point3 E(539.52898179930912, 332.45155212665065, 109.660400390625);

  cdt.insert_constraint(A, B);
  cdt.insert_constraint(C, A);
  cdt.insert_constraint(D, B);
  cdt.insert_constraint(C, E);
  return 0;
}
