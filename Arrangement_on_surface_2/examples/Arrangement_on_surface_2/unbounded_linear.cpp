#include "arr_linear.h"
#include "arr_print.h"

#include "CGAL/draw_arrangement_2.h"

int main() {
  Arrangement arr;
  auto& traits = *arr.traits();

  // Insert a n*n grid, each cell is a square of size 5
  int n = 5;
  for(int i = 0; i < n; ++i) {
    Point p1(i * 5, 0);
    Point p2(i * 5, 1);
    CGAL::insert(arr, X_monotone_curve(Line(p1, p2)));
  }

  for(int i = 0; i < n; ++i) {
    Point p1(0, i * 5);
    Point p2(1, i * 5);
    CGAL::insert(arr, X_monotone_curve(Line(p1, p2)));
  }

  Vertex_const_handle vh;
  // Generate a inner square(2*2) for all cells
  // And an inner triangle for each square
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      Point p1(i * 5 + 1, j * 5 + 1);
      Point p2(i * 5 + 4, j * 5 + 4);
      CGAL::insert(arr, X_monotone_curve(Segment(p1, Point(p2.x(), p1.y()))));
      CGAL::insert(arr, X_monotone_curve(Segment(Point(p1.x(), p2.y()), p2)));
      CGAL::insert(arr, X_monotone_curve(Segment(p1, Point(p1.x(), p2.y()))));
      CGAL::insert(arr, X_monotone_curve(Segment(Point(p2.x(), p1.y()), p2)));
      // Insert a triangle inside the square
      Point tri_p1(i * 5 + 2, j * 5 + 2);
      Point tri_p2(i * 5 + 3, j * 5 + 2);
      Point tri_p3(i * 5 + 2.5, j * 5 + 3);
      CGAL::insert(arr, X_monotone_curve(Segment(tri_p1, tri_p2)));
      CGAL::insert(arr, X_monotone_curve(Segment(tri_p2, tri_p3)));
      CGAL::insert(arr, X_monotone_curve(Segment(tri_p3, tri_p1)));
      // Connect the triangle to the square
      Point top(i * 5 + 2.5, j * 5 + 4);
      CGAL::insert(arr, X_monotone_curve(Segment(tri_p1, top)));
    }
  }

  print_arrangement_size(arr);

  CGAL::draw(arr);
}