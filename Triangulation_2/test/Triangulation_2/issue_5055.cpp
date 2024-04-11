#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>

typedef CGAL::Epick Kernel;
typedef Kernel::Point_3 Point3;

typedef CGAL::Projection_traits_xy_3<Kernel> TriangulationTraits;
typedef CGAL::Triangulation_vertex_base_2<TriangulationTraits> VertexBase;
typedef CGAL::Constrained_triangulation_face_base_2<TriangulationTraits>
    FaceBase;
typedef CGAL::Triangulation_data_structure_2<VertexBase, FaceBase>
    TriangulationData;
typedef CGAL::Constrained_Delaunay_triangulation_2<
    TriangulationTraits, TriangulationData, CGAL::Exact_predicates_tag>
    ConstrainedTriangulation;
typedef CGAL::Constrained_triangulation_plus_2<ConstrainedTriangulation> CDT;

int main() {
  CDT cdt;
  cdt.insert_constraint(
      Point3(7.22060172459952354984e+02, 5.46431441157868448499e+02,
             1.09259033203125000000e+02),
      Point3(7.21605164369421572701e+02, 5.47286947218267187054e+02,
             1.09231933593750000000e+02));
  cdt.insert_constraint(
      Point3(4.89054338074482529919e+02, 7.34795619763430352123e+01,
             1.19170654296875000000e+02),
      Point3(4.88992971641138979066e+02, 7.06742393092105203323e+01,
             1.19197509765625000000e+02));
  cdt.insert_constraint(
      Point3(6.58294921875000000000e+02, 4.94029296875000000000e+02,
             1.08803459167480468750e+02),
      Point3(6.55601638506381050320e+02, 4.92704923416975930195e+02,
             1.08774169921875000000e+02));
  cdt.insert_constraint(
      Point3(7.48077080682167661507e+02, 5.16511846264136920581e+02,
             1.09748291015625000000e+02),
      Point3(7.46966796875000000000e+02, 5.16371093750000000000e+02,
             1.09758056640625000000e+02));
  cdt.insert_constraint(
      Point3(5.59349064841414701732e+02, 4.42600717751872252848e+02,
             1.08587158203125000000e+02),
      Point3(5.56660156250000000000e+02, 4.41277343750000000000e+02,
             1.08660644531250000000e+02));
  cdt.insert_constraint(
      Point3(6.36340512288503987293e+02, 4.83223199100345823354e+02,
             1.08523681640625000000e+02),
      Point3(6.36339843750000000000e+02, 4.83224609375000000000e+02,
             1.08523681640625000000e+02));
  cdt.insert_constraint(
      Point3(5.72431981051620141443e+02, 4.62372208487347847949e+02,
             1.08266105651855468750e+02),
      Point3(5.71998046875000000000e+02, 4.61550781250000000000e+02,
             1.08266105651855468750e+02));
  cdt.insert_constraint(
      Point3(5.54615246238433428516e+02, 3.99093747119970316817e+02,
             1.08986572265625000000e+02),
      Point3(5.53910156250000000000e+02, 3.96179687500000000000e+02,
             1.09051513671875000000e+02));
  cdt.insert_constraint(
      Point3(6.29072265625000000000e+02, 4.77775390625000000000e+02,
             1.08448486328125000000e+02),
      Point3(6.29071750798135099103e+02, 4.77776482510548703431e+02,
             1.08448486328125000000e+02));
  cdt.insert_constraint(
      Point3(4.89089689771187522638e+02, 7.46450533344111875067e+01,
             1.19143554687500000000e+02),
      Point3(4.89006004028320262478e+02, 7.14945312499999943157e+01,
             1.19197509765625000000e+02));
  cdt.insert_constraint(
      Point3(4.89054338074482529919e+02, 7.34795619763430352123e+01,
             1.19170654296875000000e+02),
      Point3(4.88992951403166102864e+02, 7.06733141447368495847e+01,
             1.19197509765625000000e+02));

  return 0;
}
