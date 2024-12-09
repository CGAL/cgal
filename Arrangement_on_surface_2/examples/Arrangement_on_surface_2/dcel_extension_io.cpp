//! \file examples/Arrangement_on_surface_2/dcel_extension_io.cpp
// Using the I/O operators for arrangements with extended DCEL records.

#include <fstream>

#include <CGAL/basic.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/IO/Arr_iostream.h>
#include <CGAL/IO/Arr_text_formatter.h>

#include "arr_exact_construction_segments.h"

enum Color {BLUE, RED, WHITE};

std::ostream& operator<<(std::ostream& os, const Color& color) {
  switch (color) {
   case BLUE: os << "BLUE";  break;
   case RED: os << "RED";   break;
   case WHITE: os << "WHITE"; break;
   default: os << "ERROR!";
  }
  return os;
}

std::istream& operator>>(std::istream& is, Color& color) {
  std::string str;
  is >> str;
  if (str == "BLUE") color = BLUE;
  else if (str == "RED") color = RED;
  else if (str == "WHITE") color = WHITE;
  return is;
}

using Ext_dcel = CGAL::Arr_extended_dcel<Traits, Color, bool, int>;
using Ext_arrangement = CGAL::Arrangement_2<Traits, Ext_dcel>;
using Formatter = CGAL::Arr_extended_dcel_text_formatter<Ext_arrangement>;

int main() {
  // Construct the arrangement containing two intersecting triangles.
  Ext_arrangement arr;

  Segment s1(Point(4, 1), Point(7, 6));
  Segment s2(Point(1, 6), Point(7, 6));
  Segment s3(Point(4, 1), Point(1, 6));
  Segment s4(Point(1, 3), Point(7, 3));
  Segment s5(Point(1, 3), Point(4, 8));
  Segment s6(Point(4, 8), Point(7, 3));

  insert_non_intersecting_curve(arr, s1);
  insert_non_intersecting_curve(arr, s2);
  insert_non_intersecting_curve(arr, s3);
  insert(arr, s4);
  insert(arr, s5);
  insert(arr, s6);

  // Go over all arrangement vertices and set their colors.
  for (auto vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
    auto degree = vit->degree();
    if (degree == 0) vit->set_data(BLUE);      // Isolated vertex
    else if (degree <= 2) vit->set_data(RED);  // Vertex represents an endpoint
    else vit->set_data(WHITE);      // Vertex represents an intersection point
  }

  // Go over all arrangement edges and set their flags.
  for (auto eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
    // Check if the halfedge has the same direction as its associated
    // segment. Note that its twin always has an opposite direction.
    auto flag = (eit->source()->point() == eit->curve().source());
    eit->set_data(flag);
    eit->twin()->set_data(! flag);
  }

  // Go over all arrangement faces and print their outer boundary and indices.
  int boundary_size;
  for (auto fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
    boundary_size = 0;
    if (! fit->is_unbounded()) {
      auto curr = fit->outer_ccb();
      do ++boundary_size;
      while (++curr != fit->outer_ccb());
    }
    fit->set_data(boundary_size);
  }

  // Write the arrangement to a file.
  std::ofstream out_file("arr_ex_dcel_io.dat");
  Formatter formatter;
  CGAL::IO::write(arr, out_file, formatter);
  out_file.close();

  // Read the arrangement from the file.
  Ext_arrangement arr2;
  std::ifstream in_file("arr_ex_dcel_io.dat");

  CGAL::IO::read(arr2, in_file, formatter);
  in_file.close();

  std::cout << "The arrangement vertices:\n";
  for (auto vit = arr2.vertices_begin(); vit != arr2.vertices_end(); ++vit)
    std::cout << '(' << vit->point() << ") - " << vit->data() << std::endl;

  return 0;
}
