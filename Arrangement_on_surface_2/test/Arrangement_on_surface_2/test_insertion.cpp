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

#define N_SEGMENTS 26

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
    
  return (0);
}

