#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Combinatorial_map_operations.h>

template<typename LCC>
bool check_number_of_cells_2(LCC& lcc, unsigned int nbv, unsigned int nbe,
                             unsigned int nbf, unsigned int nbcc)
{
  std::vector<unsigned int> nbc;
  nbc=lcc.count_all_cells();

  if (nbv!=nbc[0] || nbe!=nbc[1] || nbf!=nbc[2] || nbcc!=nbc[3])
    {
      std::cout<<"Error: the number of cells is not correct. We must have "
               <<" ("<<nbv<<", "<<nbe<<", "<<nbf<<", "<<nbcc<<") and we have ("
               <<" ("<<nbc[0]<<", "<<nbc[1]<<", "<<nbc[2]<<", "<<nbc[3]<<")."
               <<std::endl;
      return false;
    }
  
  return true;
}

template<typename LCC>
bool test_LCC_2()
{
  LCC lcc;

  typedef typename LCC::Dart_handle Dart_handle;
  typedef typename LCC::Point Point;
  typedef typename LCC::Vector Vector;

  // Construction operations
  Dart_handle dh1=lcc.make_segment(Point(0,0),Point(1,0));
  Dart_handle dh2=lcc.make_segment(Point(2,0),Point(2,1));
  Dart_handle dh3=lcc.make_segment(Point(2,2),Point(3,1));
    
  lcc.template sew<0>(dh2,dh1);
  lcc.template sew<1>(dh2,dh3);
    
  Dart_handle dh5=lcc.make_triangle(Point(5,5),Point(7,5),Point(6,6));
  Dart_handle dh6=lcc.make_triangle(Point(5,4),Point(7,4),Point(6,3));
    
  lcc.template sew<2>(dh5,dh6);
    
  Dart_handle dh7=lcc.template insert_barycenter_in_cell<1>(dh1);

  Dart_handle dh8=lcc.template insert_barycenter_in_cell<2>(dh5);
    
  Dart_handle dh9=lcc.template insert_point_in_cell<1>(dh2,Point(1,0));

  Dart_handle dh10=lcc.template insert_point_in_cell<2>(dh6,Point(6,5));

  Dart_handle dh11=lcc.insert_dangling_cell_1_in_cell_2(dh8,Point(6,5.2));

  // Removal operations
  CGAL::remove_cell<LCC,1>(lcc, dh11);

  std::vector<Dart_handle> toremove;
  for ( typename LCC::template Dart_of_cell_range<0,2>::iterator
          it=lcc.template darts_of_cell<0,2>(dh10).begin(),
          itend=lcc.template darts_of_cell<0,2>(dh10).end();
          it!=itend; ++it )
    {
      Dart_handle cur = it;
      //toremove.push_back( lcc.dart_handle(it) );
      toremove.push_back( cur );
    }
  
  for ( typename std::vector<Dart_handle>::iterator
          it=toremove.begin(), itend=toremove.end(); it!=itend; ++it )
    CGAL::remove_cell<LCC,1>(lcc, *it);
  toremove.clear();
    
  CGAL::remove_cell<LCC,1>(lcc, dh9);
    
  /*    for (...) TODO with dh8
        remove_cell<LCC,1>(it); */

  CGAL::remove_cell<LCC,1>(lcc, dh7);

  lcc.template unsew<2>(dh5);

  CGAL::remove_cell<LCC,2>(lcc, dh6);
  CGAL::remove_cell<LCC,2>(lcc, dh5);
    
  lcc.template unsew<0>(dh2);
  lcc.template unsew<1>(dh2);

  CGAL::remove_cell<LCC,1>(lcc, dh1);
  CGAL::remove_cell<LCC,1>(lcc, dh2);
  CGAL::remove_cell<LCC,1>(lcc, dh3);
    
  /*    import_from_polyhedron<LCC>(lcc,ap);

        lcc.clear();
    
        import_from_plane_graph<LCC>(lcc,ais);

        lcc.clear();*/
    
  return true;
}
