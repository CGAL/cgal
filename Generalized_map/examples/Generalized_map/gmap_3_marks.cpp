#include <CGAL/Generalized_map.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Generalized_map<3> GMap_3;
typedef GMap_3::Dart_handle Dart_handle;
typedef GMap_3::size_type size_type;

int main()
{
  GMap_3 gm;

  // 1) Reserve a mark.
  size_type amark;
  try
  {
    amark = gm.get_new_mark();
  }
  catch (GMap_3::Exception_no_more_available_mark)
  {
    std::cerr<<"No more free mark, exit."<<std::endl;
    exit(-1);
  }

  // 2) Create two tetrahedra.
  Dart_handle dh1 = gm.make_combinatorial_tetrahedron();
  Dart_handle dh2 = gm.make_combinatorial_tetrahedron();

  CGAL_assertion( gm.is_valid() );
  CGAL_assertion( gm.is_volume_combinatorial_tetrahedron(dh1) );
  CGAL_assertion( gm.is_volume_combinatorial_tetrahedron(dh2) );

  // 3) 3-sew them.
  gm.sew<3>(dh1, dh2);

  // 4) Mark the darts belonging to the first tetrahedron.

  // Mark the darts belonging to the first tetrahedron.
  for (GMap_3::Dart_of_cell_range<3>::iterator
       it(gm.darts_of_cell<3>(dh1).begin()),
       itend(gm.darts_of_cell<3>(dh1).end()); it!=itend; ++it)
    gm.mark(it, amark);

  // 4) Remove the common 2-cell between the two cubes:
  //    the two tetrahedra are merged.
  gm.remove_cell<2>(dh1);

  // 5) Thanks to the mark, we know which darts come from the first tetrahedron.
  unsigned int res=0;
  for (GMap_3::Dart_range::iterator it(gm.darts().begin()),
         itend(gm.darts().end()); it!=itend; ++it)
  {
    if ( gm.is_marked(it, amark) )
      ++res;
  }

  std::cout<<"Number of darts from the first tetrahedron: "<<res<<std::endl;
  gm.free_mark(amark);

  return EXIT_SUCCESS;
}

