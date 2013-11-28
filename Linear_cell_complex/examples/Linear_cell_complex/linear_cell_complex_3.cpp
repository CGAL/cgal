#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <iostream>
#include <algorithm>

typedef CGAL::Linear_cell_complex<3> LCC_3;
typedef LCC_3::Dart_handle           Dart_handle;
typedef LCC_3::Point                 Point;
typedef LCC_3::FT                    FT;

// Functor used to display all the vertices of a given volume.
template<class LCC> 
struct Display_vol_vertices : public std::unary_function<LCC, void>
{
  Display_vol_vertices(const LCC& alcc) : 
    lcc(alcc), 
    nb_volume(0)
  {}

  void operator() (typename LCC::Dart& d) 
  { 
    std::cout<<"Volume "<<++nb_volume<<" : ";
    for (typename LCC::template One_dart_per_incident_cell_range<0,3>::
           const_iterator it=lcc.template one_dart_per_incident_cell<0,3>
           (lcc.dart_handle(d)).begin(),
           itend=lcc.template one_dart_per_incident_cell<0,3>
           (lcc.dart_handle(d)).end();
         it!=itend; ++it)
    {
      std::cout << lcc.point(it) << "; ";
    }
    std::cout<<std::endl;
  }
private:
  const LCC& lcc;
  unsigned int nb_volume;
};

int main()
{
  LCC_3 lcc;

  // Create two tetrahedra.
  Dart_handle d1 = lcc.make_tetrahedron(Point(-1, 0, 0), Point(0, 2, 0), 
                                        Point(1, 0, 0), Point(1, 1, 2));
  Dart_handle d2 = lcc.make_tetrahedron(Point(0, 2, -1),
                                        Point(-1, 0, -1),
                                        Point(1, 0, -1),
                                        Point(1, 1, -3));

  // Display all the vertices of the lcc by iterating on the 
  // Vertex_attribute container.
  CGAL::set_ascii_mode(std::cout);
  std::cout<<"Vertices: ";
  for (LCC_3::Vertex_attribute_const_range::iterator 
         v=lcc.vertex_attributes().begin(),
         vend=lcc.vertex_attributes().end(); 
       v!=vend; ++v)
    std::cout << lcc.point_of_vertex_attribute(v) << "; ";
  std::cout<<std::endl;

  // Display the vertices of each volume by iterating on darts.
  std::for_each(lcc.one_dart_per_cell<3>().begin(),
                lcc.one_dart_per_cell<3>().end(),
                Display_vol_vertices<LCC_3>(lcc));  

  // 3-Sew the 2 tetrahedra along one facet
  lcc.sew<3>(d1, d2);

  // Display the vertices of each volume by iterating on darts.
  std::for_each(lcc.one_dart_per_cell<3>().begin(),
                lcc.one_dart_per_cell<3>().end(),
                Display_vol_vertices<LCC_3>(lcc));  

  // Translate the second tetrahedra by a given vector
  LCC_3::Vector v(3,1,1);
  for (LCC_3::One_dart_per_incident_cell_range<0,3>::iterator 
         it=lcc.one_dart_per_incident_cell<0,3>(d2).begin(),
         itend=lcc.one_dart_per_incident_cell<0,3>(d2).end();
       it!=itend; ++it)
  {
    lcc.point(it)=LCC_3::Traits::Construct_translated_point_3()
      (lcc.point(it),v);
  }

  // Display the vertices of each volume by iterating on darts.
  std::for_each(lcc.one_dart_per_cell<3>().begin(),
                lcc.one_dart_per_cell<3>().end(),
                Display_vol_vertices<LCC_3>(lcc));  

  // We display the lcc characteristics.
  std::cout<<"LCC characteristics: ";
  lcc.display_characteristics(std::cout) << ", valid=" << lcc.is_valid() 
                                         << std::endl;

  return EXIT_SUCCESS;
}

