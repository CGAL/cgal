#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <iostream>
#include <algorithm>

template<typename K>
struct mypoint : public K::Point_3
{
  typedef typename K::Point_3 Base;
  
  mypoint() : mtype('a')
  {}

  mypoint(const Base& apoint) : Base(apoint),
                                mtype('b')
  {}
  
  mypoint(const mypoint& apoint) : Base(apoint),
                                   mtype('c')
  {}

  mypoint(typename K::FT a, typename K::FT b, typename K::FT c) : Base(a,b,c),
                                                                  mtype('d')
  {}

  char type() const
  { return mtype; }
  
private:
  char mtype;
};

template<typename K>
struct mytraits: public CGAL::Linear_cell_complex_traits<3, K>
{
  typedef mypoint<K> Point;
};

typedef mytraits<CGAL::Exact_predicates_inexact_constructions_kernel> Traits;

typedef CGAL::Linear_cell_complex<3,3, Traits> LCC_3;
typedef LCC_3::Dart_handle                     Dart_handle;
typedef LCC_3::Point                           Point;
typedef LCC_3::FT                              FT;

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

  // 3-Sew the two cubes along one facet.
  lcc.sew<3>(lcc.beta(d1, 1, 1, 2), lcc.beta(d2, 2));

  // Barycentric triangulation of the facet between the two cubes.
  lcc.insert_barycenter_in_cell<2>(lcc.beta(d2, 2));

  // Display all the vertices of the map.
  for (LCC_3::Vertex_attribute_range::iterator 
         it=lcc.vertex_attributes().begin(),
         itend=lcc.vertex_attributes().end(); 
       it!=itend; ++it)
  {
    std::cout<<"point: "<<lcc.point_of_vertex_attribute(it)
             <<", "<<"type: "<<lcc.point_of_vertex_attribute(it).type()
             <<std::endl;
  }

  return EXIT_SUCCESS;
}
