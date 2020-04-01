// Testing the insertion functions, covering all possible scenarios.

#include <CGAL/Quotient.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>

typedef CGAL::Quotient<int>                           Number_type;
typedef CGAL::Simple_cartesian<Number_type>           Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>            Traits_2;
typedef Traits_2::Point_2                             Point_2;
typedef Traits_2::X_monotone_curve_2                  Segment_2;
typedef CGAL::Arrangement_2<Traits_2>                 Arrangement_2;
typedef Arrangement_2::Halfedge_handle                Halfedge_handle;
typedef Arrangement_2::Vertex_handle                  Vertex_handle;

#define N_SEGMENTS 26

//test from a bug report when inserting a segment closing a hole
bool test_insert_at_vertices_1(){
  Arrangement_2 arr;

  Vertex_handle v_x0y3 = arr.insert_in_face_interior(Point_2(Number_type(0), Number_type(3)), arr.unbounded_face());
  Vertex_handle v_x1y6 = arr.insert_in_face_interior(Point_2(Number_type(1), Number_type(6)), arr.unbounded_face());
  Vertex_handle v_x1y3 = arr.insert_in_face_interior(Point_2(Number_type(1), Number_type(3)), arr.unbounded_face());
  Vertex_handle v_x2y3 = arr.insert_in_face_interior(Point_2(Number_type(2), Number_type(3)), arr.unbounded_face());
  Vertex_handle v_x3y6 = arr.insert_in_face_interior(Point_2(Number_type(3), Number_type(6)), arr.unbounded_face());
  Vertex_handle v_x3y0 = arr.insert_in_face_interior(Point_2(Number_type(3), Number_type(0)), arr.unbounded_face());

  arr.insert_at_vertices(Segment_2(v_x0y3->point(), v_x1y6->point()), v_x0y3, v_x1y6);
  arr.insert_at_vertices(Segment_2(v_x0y3->point(), v_x1y3->point()), v_x0y3, v_x1y3);
  arr.insert_at_vertices(Segment_2(v_x1y3->point(), v_x2y3->point()), v_x1y3, v_x2y3);
  arr.insert_at_vertices(Segment_2(v_x0y3->point(), v_x3y6->point()), v_x0y3, v_x3y6);
  arr.insert_at_vertices(Segment_2(v_x0y3->point(), v_x3y0->point()), v_x0y3, v_x3y0);

  Halfedge_handle he = arr.insert_at_vertices(Segment_2(v_x3y6->point(), v_x3y0->point()), v_x3y6, v_x3y0);


  if (he->face() != arr.unbounded_face())
  {
    std::cerr << "Error: he->face() must be the unbounded face!" << std::endl;
    return false;
  }

  return is_valid(arr);
}

bool test_insert_at_vertices_2(){
  Arrangement_2 arr;

  Vertex_handle v_x0y3 = arr.insert_in_face_interior(Kernel::Point_2(Kernel::FT(0), Kernel::FT(3)), arr.unbounded_face());
  Vertex_handle v_x1y6 = arr.insert_in_face_interior(Kernel::Point_2(Kernel::FT(1), Kernel::FT(6)), arr.unbounded_face());
  Vertex_handle v_x1y3 = arr.insert_in_face_interior(Kernel::Point_2(Kernel::FT(1), Kernel::FT(3)), arr.unbounded_face());
  Vertex_handle v_x2y3 = arr.insert_in_face_interior(Kernel::Point_2(Kernel::FT(2), Kernel::FT(3)), arr.unbounded_face());
  Vertex_handle v_x3y6 = arr.insert_in_face_interior(Kernel::Point_2(Kernel::FT(3), Kernel::FT(6)), arr.unbounded_face());
  Vertex_handle v_x3y0 = arr.insert_in_face_interior(Kernel::Point_2(Kernel::FT(3), Kernel::FT(0)), arr.unbounded_face());

  arr.insert_at_vertices(Segment_2(v_x0y3->point(), v_x1y6->point()), v_x0y3, v_x1y6);
  arr.insert_at_vertices(Segment_2(v_x0y3->point(), v_x1y3->point()), v_x0y3, v_x1y3);
  arr.insert_at_vertices(Segment_2(v_x1y3->point(), v_x2y3->point()), v_x1y3, v_x2y3);
  arr.insert_at_vertices(Segment_2(v_x0y3->point(), v_x3y6->point()), v_x0y3, v_x3y6);
  arr.insert_at_vertices(Segment_2(v_x3y6->point(), v_x3y0->point()), v_x3y6, v_x3y0);

  Halfedge_handle he = arr.insert_at_vertices(Segment_2(v_x3y0->point(), v_x0y3->point()), v_x3y0, v_x0y3);

  if (he->face() != arr.unbounded_face())
  {
    std::cerr << "Error: he->face() must be the unbounded face!" << std::endl;
    return false;
  }

  return is_valid(arr);
}

bool test_insert_at_vertices(){
  return test_insert_at_vertices_1() && test_insert_at_vertices_2();
}

int main ()
{
  Arrangement_2   arr;
  Segment_2       segs[N_SEGMENTS];
  bool            valid;
  int             k;

  segs[0] = Segment_2 (Point_2 (5, 9), Point_2 (5, 11));
  segs[1] = Segment_2 (Point_2 (5, 9), Point_2 (7, 9));
  segs[2] = Segment_2 (Point_2 (5, 11), Point_2 (7, 9));
  segs[3] = Segment_2 (Point_2 (7, 6), Point_2 (9, 7));
  segs[4] = Segment_2 (Point_2 (9, 7), Point_2 (7, 9));
  segs[5] = Segment_2 (Point_2 (5, 9), Point_2 (7, 6));
  segs[6] = Segment_2 (Point_2 (10, 11), Point_2 (8, 13));
  segs[7] = Segment_2 (Point_2 (8, 13), Point_2 (11, 12));
  segs[8] = Segment_2 (Point_2 (10, 11), Point_2 (11, 12));
  segs[9] = Segment_2 (Point_2 (1, 20), Point_2 (5, 1));
  segs[10] = Segment_2 (Point_2 (5, 1), Point_2 (12, 6));
  segs[11] = Segment_2 (Point_2 (1, 20), Point_2 (12, 6));
  segs[12] = Segment_2 (Point_2 (13, 13), Point_2 (13, 15));
  segs[13] = Segment_2 (Point_2 (13, 15), Point_2 (15, 12));
  segs[14] = Segment_2 (Point_2 (13, 13), Point_2 (15, 12));
  segs[15] = Segment_2 (Point_2 (11, 12), Point_2 (13, 13));
  segs[16] = Segment_2 (Point_2 (1, 20), Point_2 (11, 17));
  segs[17] = Segment_2 (Point_2 (11, 17), Point_2 (20, 13));
  segs[18] = Segment_2 (Point_2 (12, 6), Point_2 (20, 13));
  segs[19] = Segment_2 (Point_2 (15, 12), Point_2 (20, 13));
  segs[20] = Segment_2 (Point_2 (5, 1), Point_2 (20, 1));
  segs[21] = Segment_2 (Point_2 (20, 1), Point_2 (20, 13));
  segs[22] = Segment_2 (Point_2 (15, 5), Point_2 (17, 3));
  segs[23] = Segment_2 (Point_2 (13, 3), Point_2 (17, 3));
  segs[24] = Segment_2 (Point_2 (12, 6), Point_2 (13, 3));
  segs[25] = Segment_2 (Point_2 (17, 3), Point_2 (20, 1));

  for (k = 0; k < N_SEGMENTS; k++)
  {
    insert_non_intersecting_curve (arr, segs[k]);
    valid = arr.is_valid();
    std::cout << "  Inserted " << k+1 << " segment(s), arrangement is "
              << (valid ? "valid." : "NOT valid!") << std::endl;

    if (! valid)
      return (1);
  }

  std::cout << "Arrangement size:"
            << "   V = " << arr.number_of_vertices()
            << ",  E = " << arr.number_of_edges()
            << ",  F = " << arr.number_of_faces() << std::endl;

  // Check the validity more thoroughly.
  valid = is_valid(arr);
  std::cout << "Arrangement is "
            << (valid ? "valid." : "NOT valid!") << std::endl;
  if (!valid) return 1;

  std::cout << "Test insert_at_vertices ";
  valid=test_insert_at_vertices();
  std::cout << ( valid ? "valid." : "NOT valid!") << std::endl;

  return valid?0:1;
}

