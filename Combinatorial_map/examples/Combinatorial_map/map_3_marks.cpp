#include <CGAL/Combinatorial_map.h>
#include <CGAL/Combinatorial_map_constructors.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <iostream>
#include <cstdlib>

typedef CGAL::Combinatorial_map<3> CMap_3;
typedef CMap_3::Dart_handle Dart_handle;

int main()
{
  CMap_3 cm;
  
  // 1) Reserve a mark.
  int mark = cm.get_new_mark();
  if ( mark==-1 )
    {
      std::cerr<<"No more free mark, exit."<<std::endl;
      exit(-1);
    }
  
  // 2) Create two tetrahedra.
  Dart_handle dh1 = CGAL::make_combinatorial_tetrahedron(cm);  
  Dart_handle dh2 = CGAL::make_combinatorial_tetrahedron(cm);

  // 3) 3-sew them.
  cm.sew<3>(dh1, dh2);
  
  // 4) Mark the darts belonging to the first tetrahedron.
  for  (CMap_3::Dart_of_cell_range<3>::iterator 
          it(cm.darts_of_cell<3>(dh1).begin()),
          itend(cm.darts_of_cell<3>(dh1).end()); it!=itend; ++it)
    cm.mark(it, mark);

  // 4) Remove the common 2-cell between the two cubes:
  // the two tetrahedra are merged.
  CGAL::remove_cell<CMap_3, 2>(cm, dh1);

  // 5) Thanks to the mark, we know which darts come from the first tetrahedron.
  unsigned int res=0;
  for (CMap_3::Dart_range::iterator it(cm.darts().begin()),
	 itend(cm.darts().end()); it!=itend; ++it)
    {
      if ( cm.is_marked(it, mark) )
	++res;
    }
  
  std::cout<<"Number of darts from the first tetrahedron: "<<res<<std::endl;
  cm.free_mark(mark);

  return EXIT_SUCCESS;
}

