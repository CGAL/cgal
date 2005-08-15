// file: examples/Arrangement_2/example20.C

#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_face_extended_dcel.h>

typedef int                                             Number_type;
typedef CGAL::Cartesian<Number_type>                    Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>              Traits_2;
typedef Traits_2::Point_2                               Point_2;
typedef Traits_2::X_monotone_curve_2                    Segment_2;
typedef CGAL::Arr_face_extended_dcel<int, Traits_2>     Dcel;
typedef CGAL::Arrangement_2<Traits_2, Dcel>             Arrangement_2;
typedef Arrangement_2::Vertex_handle                    Vertex_handle;
typedef Arrangement_2::Halfedge_handle                  Halfedge_handle;

int main()
{
  Arrangement_2   arr;

  Segment_2 s1(Point_2(1, 3), Point_2(3, 5));
  Segment_2 s2(Point_2(3, 5), Point_2(5, 3));
  Segment_2 s3(Point_2(5, 3), Point_2(3, 1));
  Segment_2 s4(Point_2(3, 1), Point_2(1, 3));
  Segment_2 s5(Point_2(1, 3), Point_2(5, 3));

  Halfedge_handle e1 = arr.insert_in_face_interior(s1, arr.unbounded_face());
  Vertex_handle   v1 = e1->source();
  Vertex_handle   v2 = e1->target();
  Halfedge_handle e2 = arr.insert_from_left_vertex(s2, v2);
  Vertex_handle   v3 = e2->target();
  Halfedge_handle e3 = arr.insert_from_right_vertex(s3, v3);
  Vertex_handle   v4 = e3->target();
  Halfedge_handle e4 = arr.insert_at_vertices(s4, v4, v1);
  Halfedge_handle e5 = arr.insert_at_vertices(s5, v1, v3);

  int i = 0;
  Arrangement_2::Face_iterator fit;
  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
    fit->set_data(i++);

  // Print all neighboring faces to all faces
  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    if(fit->is_unbounded()) continue;          // skip the unbounded face
    std::cout << fit->get_data() << ":" << std::endl;
    Arrangement_2::Ccb_halfedge_circulator hec = fit->outer_ccb();
    Arrangement_2::Ccb_halfedge_circulator hec_begin = hec;

    Arrangement_2::Face_handle neighbor_face = hec->twin()->face();
    std::cout << neighbor_face->get_data();
    for (++hec; hec != hec_begin; ++hec) { 
      neighbor_face = hec->twin()->face();
      std::cout << ", " << neighbor_face->get_data();
    }
    std::cout << std::endl;
  }
  return 0;
}  
