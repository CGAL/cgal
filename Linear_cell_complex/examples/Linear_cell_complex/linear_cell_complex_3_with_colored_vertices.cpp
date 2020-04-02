#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <iostream>
#include <algorithm>

struct Average_functor
{
  template<class CellAttribute>
  void operator()(CellAttribute& ca1, CellAttribute& ca2)
  { ca1.info()=(ca1.info()+ ca2.info())/2; }
};

struct Myitem
{
  template<class Refs>
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute_with_point< Refs, int, CGAL::Tag_true,
                                             Average_functor >
    Vertex_attribute;

    typedef std::tuple<Vertex_attribute> Attributes;
  };
};

typedef CGAL::Linear_cell_complex_traits
<3, CGAL::Exact_predicates_inexact_constructions_kernel> Traits;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3,Traits,Myitem> LCC_3;
typedef LCC_3::Dart_handle                               Dart_handle;
typedef LCC_3::Point                                     Point;
typedef LCC_3::FT                                        FT;

Dart_handle make_iso_cuboid(LCC_3& lcc, const Point& basepoint, FT lg)
{
  return lcc.make_hexahedron(basepoint,
                             Traits::Construct_translated_point()
                             (basepoint,Traits::Vector(lg,0,0)),
                             Traits::Construct_translated_point()
                             (basepoint,Traits::Vector(lg,lg,0)),
                             Traits::Construct_translated_point()
                             (basepoint,Traits::Vector(0,lg,0)),
                             Traits::Construct_translated_point()
                             (basepoint,Traits::Vector(0,lg,lg)),
                             Traits::Construct_translated_point()
                             (basepoint,Traits::Vector(0,0,lg)),
                             Traits::Construct_translated_point()
                             (basepoint,Traits::Vector(lg,0,lg)),
                             Traits::Construct_translated_point()
                             (basepoint,Traits::Vector(lg,lg,lg)));
}

int main()
{
  LCC_3 lcc;

  // Create two iso_cuboids.
  Dart_handle d1 = make_iso_cuboid(lcc, Point(-2, 0, 0), 1);
  Dart_handle d2 = make_iso_cuboid(lcc, Point(0, 0, 0), 1);

  // Set the "color" of all vertices of the first cube to 1.
  for (LCC_3::One_dart_per_incident_cell_range<0, 3>::iterator
         it=lcc.one_dart_per_incident_cell<0,3>(d1).begin(),
         itend=lcc.one_dart_per_incident_cell<0,3>(d1).end(); it!=itend; ++it)
  { lcc.info<0>(it)=1; }

  // Set the "color" of all vertices of the second cube to 19.
  for (LCC_3::One_dart_per_incident_cell_range<0, 3>::iterator it=
         lcc.one_dart_per_incident_cell<0,3>(d2).begin(),
         itend=lcc.one_dart_per_incident_cell<0,3>(d2).end(); it!=itend; ++it)
  { lcc.info<0>(it)=19; }

  // 3-Sew the two cubes along one facet.
  lcc.sew<3>(lcc.beta(d1, 1, 1, 2), lcc.beta(d2, 2));

  // Barycentric triangulation of the facet between the two cubes.
  Dart_handle d3=lcc.insert_barycenter_in_cell<2>(lcc.beta(d2, 2));

  // Set the color of the new vertex to 5.
  lcc.info<0>(d3)=5;

  // Display all the vertices of the map.
  for (LCC_3::Vertex_attribute_range::iterator
         it=lcc.vertex_attributes().begin(),
         itend=lcc.vertex_attributes().end();
       it!=itend; ++it)
  {
    std::cout<<"point: "<<lcc.point_of_vertex_attribute(it)<<", "<<"color: "
             <<lcc.info_of_attribute<0>(it)<<std::endl;
  }

  return EXIT_SUCCESS;
}
