#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <iostream>
#include <algorithm>

typedef CGAL::Combinatorial_map<3> CMap_3;
typedef CMap_3::Dart_handle        Dart_handle;

// Functor used to display all the vertices of a given volume
template<class CMap, unsigned int i>
struct Display_vertices_of_cell : public CGAL::cpp98::unary_function<CMap, void>
{
  Display_vertices_of_cell(const CMap& acmap) :
    cmap(acmap),
    nb_cell(0)
  {}

  void operator() (const typename CMap::Dart& d)
  {
    std::cout<<i<<"-cell "<<++nb_cell<<" : ";
    for (typename CMap::template One_dart_per_incident_cell_range<0,i>::
           const_iterator it=cmap.template one_dart_per_incident_cell<0,i>
           (cmap.dart_handle(d)).begin(),
           itend=cmap.template one_dart_per_incident_cell<0,i>
           (cmap.dart_handle(d)).end();
         it!=itend; ++it)
    {
      std::cout << cmap.darts().index(it) << "; ";
    }
    std::cout<<std::endl;
  }

  void operator() (const typename CMap::Dart* ptr)
  { operator() (*ptr); }

private:
  const CMap& cmap;
  unsigned int nb_cell;
};

// Functor used to remove a face
template<class CMap>
struct Remove_face : public CGAL::cpp98::unary_function<CMap, void>
{
  Remove_face(CMap& acmap) : cmap(acmap)
  {}

  void operator() (typename CMap::Dart* d)
  {
    cmap.template remove_cell<2>(cmap.dart_handle(*d));
    std::cout<<"CMap characteristics: ";
    cmap.display_characteristics(std::cout) << ", valid=" << cmap.is_valid()
                                            << std::endl;
  }

private:
  CMap& cmap;
};


// Functor allowing to transform a variable into its address.
template<typename T>
struct Take_address : public CGAL::cpp98::unary_function<T, T*>
{
  T* operator() (T& t) const
  { return &t; }
};

int main()
{
  CMap_3 cmap;

  // Create two tetrahedra.
  Dart_handle d1 = cmap.make_combinatorial_tetrahedron();
  Dart_handle d2 = cmap.make_combinatorial_tetrahedron();

  // Display the vertices of each volume by iterating on darts.
  std::cout<<"********Volumes********"<<std::endl;
  std::for_each(cmap.one_dart_per_cell<3>().begin(),
                cmap.one_dart_per_cell<3>().end(),
                Display_vertices_of_cell<CMap_3,3>(cmap));

  // 3-Sew the 2 tetrahedra along one facet
  cmap.sew<3>(d1, d2);

  // Display the vertices of each face by iterating on darts.
  std::cout<<"********Faces********"<<std::endl;
  std::for_each(cmap.one_dart_per_cell<2>().begin(),
                cmap.one_dart_per_cell<2>().end(),
                Display_vertices_of_cell<CMap_3,2>(cmap));

    // We display the map characteristics.
  std::cout<<"CMap characteristics: ";
  cmap.display_characteristics(std::cout) << ", valid=" << cmap.is_valid()
                                          << std::endl << std::endl;

  std::vector<CMap_3::Dart*> toremove;

  // Copy in vector toremove one dart per face
  std::copy(boost::transform_iterator<Take_address<CMap_3::Dart>,
                                      CMap_3::One_dart_per_cell_range<2>::iterator>
            (cmap.one_dart_per_cell<2>().begin(),
             Take_address<CMap_3::Dart>()),
            boost::transform_iterator<Take_address<CMap_3::Dart>,
                                      CMap_3::One_dart_per_cell_range<2>::iterator>
            (cmap.one_dart_per_cell<2>().end(),
             Take_address<CMap_3::Dart>()),
            back_inserter(toremove));

  // Remove each face sequentially.
  std::for_each(toremove.begin(), toremove.end(), Remove_face<CMap_3>(cmap));

  return EXIT_SUCCESS;
}
