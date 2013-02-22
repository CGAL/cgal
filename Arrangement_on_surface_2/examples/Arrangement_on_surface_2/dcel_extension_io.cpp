//! \file examples/Arrangement_on_surface_2/dcel_extension_io.cpp
// Using the I/O operators for arrangements with extended DCEL records.

#include "arr_rational_nt.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/IO/Arr_text_formatter.h>
#include <CGAL/IO/Arr_iostream.h>
#include <fstream>

enum Color {BLUE, RED, WHITE};

std::ostream& operator<< (std::ostream& os, const Color& color)
{
  switch (color)
  {
  case BLUE:  os << "BLUE";  break;
  case RED:   os << "RED";   break;
  case WHITE: os << "WHITE"; break;
  default: os << "ERROR!";
  }
  return (os);
}

std::istream& operator>> (std::istream& is, Color& color)
{
  std::string   str;
  is >> str;

  if (str == "BLUE")
    color = BLUE;
  else if (str == "RED")
    color = RED;
  else if (str == "WHITE")
    color = WHITE;

  return (is);
}

typedef CGAL::Cartesian<Number_type>                      Kernel;
typedef CGAL::Arr_segment_traits_2<Kernel>                Traits_2;
typedef Traits_2::Point_2                                 Point_2;
typedef Traits_2::X_monotone_curve_2                      Segment_2;
typedef CGAL::Arr_extended_dcel<Traits_2,
                                Color, bool, int>         Dcel;
typedef CGAL::Arrangement_2<Traits_2, Dcel>               Arrangement_2;
typedef CGAL::Arr_extended_dcel_text_formatter<Arrangement_2>  Formatter;

int main ()
{
  // Construct the arrangement containing two intersecting triangles.
  Arrangement_2          arr;

  Segment_2      s1 (Point_2(4, 1), Point_2(7, 6));
  Segment_2      s2 (Point_2(1, 6), Point_2(7, 6));
  Segment_2      s3 (Point_2(4, 1), Point_2(1, 6));
  Segment_2      s4 (Point_2(1, 3), Point_2(7, 3));
  Segment_2      s5 (Point_2(1, 3), Point_2(4, 8));
  Segment_2      s6 (Point_2(4, 8), Point_2(7, 3));

  insert_non_intersecting_curve (arr, s1);
  insert_non_intersecting_curve (arr, s2);
  insert_non_intersecting_curve (arr, s3);
  insert (arr, s4);
  insert (arr, s5);
  insert (arr, s6);

  // Go over all arrangement vertices and set their colors.
  Arrangement_2::Vertex_iterator            vit;
  std::size_t                               degree;

  for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit)
  {
    degree = vit->degree();
    if (degree == 0)
      vit->set_data (BLUE);       // Isolated vertex.
    else if (degree <= 2)
      vit->set_data (RED);        // Vertex represents an endpoint.
    else
      vit->set_data (WHITE);      // Vertex represents an intersection point.
  }

  // Go over all arrangement edges and set their flags.
  Arrangement_2::Edge_iterator              eit;
  bool                                      flag;

  for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit)
  {
    // Check if the halfegde has the same diretion as its associated
    // segment. Note that its twin always has an opposite direction.
    flag = (eit->source()->point() == eit->curve().source());
    eit->set_data (flag);
    eit->twin()->set_data (!flag);
  }

  // Go over all arrangement faces and print their outer boundary and indices.
  Arrangement_2::Face_iterator              fit;
  Arrangement_2::Ccb_halfedge_circulator    curr;
  int                                       boundary_size;

  for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
  {
    boundary_size = 0;
    if (! fit->is_unbounded())
    {
      curr = fit->outer_ccb();
      do
      {
        ++boundary_size;
        ++curr;
      } while (curr != fit->outer_ccb());
    }
    fit->set_data (boundary_size);
  }

  // Write the arrangement to a file.
  std::ofstream    out_file ("arr_ex_dcel_io.dat");
  Formatter        formatter;

  write (arr, out_file, formatter);
  out_file.close();

  // Read the arrangement from the file.
  Arrangement_2    arr2;
  std::ifstream    in_file ("arr_ex_dcel_io.dat");

  read (arr2, in_file, formatter);
  in_file.close();

  std::cout << "The arrangement vertices: " << std::endl;
  for (vit = arr2.vertices_begin(); vit != arr2.vertices_end(); ++vit)
    std::cout << '(' << vit->point() << ") - " << vit->data() << std::endl;

  return (0);
}
