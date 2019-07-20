#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/natural_neighbor_coordinates_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

#include <iostream>
#include <iterator>
#include <utility>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel          K;

typedef K::FT                                                        Coord_type;
typedef CGAL::Triangulation_vertex_base_with_info_2<Coord_type, K>   Vb;
typedef CGAL::Triangulation_data_structure_2<Vb>                     Tds;
typedef CGAL::Delaunay_triangulation_2<K, Tds>                       Delaunay_triangulation;

// The functor 'Identity' matches anything to itself
typedef Delaunay_triangulation::Vertex_handle                        Vertex_handle;
typedef CGAL::Identity<std::pair<Vertex_handle, Coord_type> >        Identity;

// Resulting points-coordinates pairs are here stored in an object of this type:
typedef std::vector<std::pair<Vertex_handle, Coord_type> >           Point_coordinate_vector;

int main()
{
  Delaunay_triangulation dt;

  for(int y=0; y<3; ++y)
    for(int x=0; x<3; ++x)
      dt.insert(K::Point_2(x, y));

  // coordinates computation
  K::Point_2 p(1.2, 0.7); // query point
  Point_coordinate_vector coords;

  // The functor Identity is passed to the method
  CGAL::Triple<std::back_insert_iterator<Point_coordinate_vector>, K::FT, bool> result =
      CGAL::natural_neighbor_coordinates_2(dt, p, std::back_inserter(coords), Identity());

  if(!result.third)
  {
    std::cout << "The coordinate computation was not successful." << std::endl;
    std::cout << "The point (" << p << ") lies outside the convex hull." << std::endl;
  }

  // Assign the coordinates to the vertices
  std::cout << "Coordinates for point: (" << p << ") are the following: " << std::endl;
  for(std::size_t i=0; i<coords.size(); ++i)
  {
    Vertex_handle vh = coords[i].first;
    vh->info() = coords[i].second;
    std::cout << "  Vertex: (" << vh->point() << ") coeff: " << vh->info() << std::endl;
  }

  return EXIT_SUCCESS;
}
