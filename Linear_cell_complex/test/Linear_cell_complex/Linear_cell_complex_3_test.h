#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Combinatorial_map_operations.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <fstream>

template<typename LCC>
bool check_number_of_cells_3(LCC& lcc, unsigned int nbv, unsigned int nbe,
                             unsigned int nbf, unsigned int nbvol,
                             unsigned int nbcc)
{
  if ( !lcc.is_valid() )
  {
    std::cout<<"ERROR: the lcc is not valid."<<std::endl;
    assert(false);
    return false;
  }
  
  std::vector<unsigned int> nbc;
  nbc=lcc.count_all_cells();

  if (nbv!=nbc[0] || nbe!=nbc[1] || nbf!=nbc[2] || nbvol!=nbc[3] ||
      nbcc!=nbc[4])
  {
    std::cout<<"ERROR: the number of cells is not correct. We must have "
             <<" ("<<nbv<<", "<<nbe<<", "<<nbf<<", "<<nbvol<<", "<<nbcc
             <<") and we have"<<" ("<<nbc[0]<<", "<<nbc[1]<<", "<<nbc[2]<<", "
             <<nbc[3]<<", "<<nbc[4]<<")."
             <<std::endl;
    assert(false);
    return false;
  }

  if ( nbv!=lcc.number_of_vertex_attributes() )
  {
    std::cout<<"ERROR: the number of vertices ("<<nbv<<") is different than "
             <<"the number of vertex attributes ("
             <<lcc.number_of_vertex_attributes()<<")"<<std::endl;

    assert(false);
    return false;
  }
  
  return true;
}

template<typename LCC>
bool test_LCC_3()
{
  LCC lcc;

  typedef typename LCC::Dart_handle Dart_handle;
  typedef typename LCC::Point Point;
  typedef typename LCC::Vector Vector;
  
  // Construction operations
  Dart_handle dh1=lcc.make_segment(Point(0,0,0),Point(1,0,0));
  Dart_handle dh2=lcc.make_segment(Point(2,0,0),Point(2,1,0));
  Dart_handle dh3=lcc.make_segment(Point(2,2,0),Point(3,1,0));
  if ( !check_number_of_cells_3(lcc, 6, 3, 6, 3, 3) )
    return false;
  
  lcc.template sew<0>(dh2,dh1);
  lcc.template sew<1>(dh2,dh3);
  if ( !check_number_of_cells_3(lcc, 4, 3, 4, 1, 1) )
    return false;

  Dart_handle dh5=lcc.make_triangle(Point(5,5,3),Point(7,5,3),Point(6,6,3));
  Dart_handle dh6=lcc.make_triangle(Point(5,4,3),Point(7,4,3),Point(6,3,3));    
  if ( !check_number_of_cells_3(lcc, 10, 9, 6, 3, 3) )
    return false;

  lcc.template sew<2>(dh5,dh6);
  if ( !check_number_of_cells_3(lcc, 8, 8, 6, 2, 2) )
    return false;

  Dart_handle dh7=lcc.template insert_barycenter_in_cell<1>(dh1);
  if ( !check_number_of_cells_3(lcc, 9, 9, 6, 2, 2) )
    return false;

  Dart_handle dh8=lcc.template insert_barycenter_in_cell<2>(dh5);
  if ( !check_number_of_cells_3(lcc, 10, 12, 8, 2, 2) )
    return false;

  Dart_handle dh9=lcc.template insert_point_in_cell<1>(dh2,Point(1,0,3));
  if ( !check_number_of_cells_3(lcc, 11, 13, 8, 2, 2) )
    return false;

  Dart_handle dh10=lcc.template insert_point_in_cell<2>(dh6,Point(6,5,3));
  if ( !check_number_of_cells_3(lcc, 12, 16, 10, 2, 2) )
    return false;

  Dart_handle dh11=lcc.insert_dangling_cell_1_in_cell_2(dh8,Point(6,5.2,3));
  if ( !check_number_of_cells_3(lcc, 13, 17, 10, 2, 2) )
    return false;

  Dart_handle dh12 = lcc.make_tetrahedron(Point(-1, 0, 0),Point(0, 2, 0),
                                          Point(1, 0, 0),Point(1, 1, 2));
  Dart_handle dh13 = lcc.make_tetrahedron(Point(0, 2, -1),Point(-1, 0, -1),
                                          Point(1, 0, -1),Point(1, 1, -3));
  if ( !check_number_of_cells_3(lcc, 21, 29, 18, 4, 4) )
    return false;

  lcc.template sew<3>(dh12, dh13);
  if ( !check_number_of_cells_3(lcc, 18, 26, 17, 4, 3) )
    return false;

  Dart_handle dh14=lcc.template insert_barycenter_in_cell<2>(dh12);
  if ( !check_number_of_cells_3(lcc, 19, 29, 19, 4, 3) )
    return false;

  Dart_handle dh15=lcc.template insert_barycenter_in_cell<1>(dh14);
  if ( !check_number_of_cells_3(lcc, 20, 30, 19, 4, 3) )
    return false;

  // Removal operations
  CGAL::remove_cell<LCC,0>(lcc, dh15);
  if ( !check_number_of_cells_3(lcc, 19, 29, 19, 4, 3) )
    return false;

  CGAL::remove_cell<LCC,1>(lcc, dh14->beta(2)->beta(1));
  CGAL::remove_cell<LCC,1>(lcc, dh14->beta(0));
  CGAL::remove_cell<LCC,1>(lcc, dh14);
  if ( !check_number_of_cells_3(lcc, 18, 26, 17, 4, 3) )
    return false;
  
  lcc.template unsew<3>(dh12);
  if ( !check_number_of_cells_3(lcc, 21, 29, 18, 4, 4) )
    return false;

  CGAL::remove_cell<LCC,3>(lcc, dh13);
  CGAL::remove_cell<LCC,3>(lcc, dh12);
  if ( !check_number_of_cells_3(lcc, 13, 17, 10, 2, 2) )
    return false;
  
  CGAL::remove_cell<LCC,1>(lcc, dh11);
  if ( !check_number_of_cells_3(lcc, 12, 16, 10, 2, 2) )
    return false;

  std::vector<Dart_handle> toremove;
  for ( typename LCC::template Dart_of_cell_range<0,2>::iterator
          it=lcc.template darts_of_cell<0,2>(dh10).begin(),
          itend=lcc.template darts_of_cell<0,2>(dh10).end();
        it!=itend; ++it )
    toremove.push_back( it );
  
  for ( typename std::vector<Dart_handle>::iterator
          it=toremove.begin(), itend=toremove.end(); it!=itend; ++it )
    CGAL::remove_cell<LCC,1>(lcc, *it);
  toremove.clear();
  if ( !check_number_of_cells_3(lcc, 11, 13, 8, 2, 2) )
    return false;
  
  CGAL::remove_cell<LCC,0>(lcc, dh9);
  if ( !check_number_of_cells_3(lcc, 10, 12, 8, 2, 2) )
    return false;

  for ( typename LCC::template Dart_of_cell_range<0,2>::iterator
          it=lcc.template darts_of_cell<0,2>(dh8).begin(),
          itend=lcc.template darts_of_cell<0,2>(dh8).end();
        it!=itend; ++it )
    toremove.push_back( it );
  
  for ( typename std::vector<Dart_handle>::iterator
          it=toremove.begin(), itend=toremove.end(); it!=itend; ++it )
    CGAL::remove_cell<LCC,1>(lcc, *it);
  toremove.clear();
  if ( !check_number_of_cells_3(lcc, 9, 9, 6, 2, 2) )
    return false;
  
  CGAL::remove_cell<LCC,0>(lcc, dh7);
  if ( !check_number_of_cells_3(lcc, 8, 8, 6, 2, 2) )
    return false;

  lcc.template unsew<2>(dh5);
  if ( !check_number_of_cells_3(lcc, 10, 9, 6, 3, 3) )
    return false;

  CGAL::remove_cell<LCC,2>(lcc, dh6);
  CGAL::remove_cell<LCC,2>(lcc, dh5);
  if ( !check_number_of_cells_3(lcc, 4, 3, 4, 1, 1) )
    return false;

  lcc.template unsew<1>(dh2);
  if ( !check_number_of_cells_3(lcc, 5, 3, 5, 2, 2) )
    return false;

  lcc.template unsew<0>(dh2); 
  if ( !check_number_of_cells_3(lcc, 6, 3, 6, 3, 3) )
    return false;

  CGAL::remove_cell<LCC,1>(lcc, dh1);
  CGAL::remove_cell<LCC,1>(lcc, dh2);
  CGAL::remove_cell<LCC,1>(lcc, dh3);  
  if ( !check_number_of_cells_3(lcc, 0, 0, 0, 0, 0) )
    return false;


  {
    CGAL::Polyhedron_3<typename LCC::Traits> P;
    std::ifstream in("data/armadillo.off");
    if ( in.fail() )
    {
      std::cout<<"Error: impossible to open 'data/armadillo.off'"<<std::endl;
      return false;
    }
    in >> P;
    CGAL::import_from_polyhedron_3<LCC>(lcc,P);
    if ( !check_number_of_cells_3(lcc, 26002, 78000, 52000, 1, 1) )
      return false;
    lcc.clear();
  }

  {
    CGAL::Triangulation_3<typename LCC::Traits> T;
    std::ifstream in("data/points.txt");
    if ( in.fail() )
    {
      std::cout<<"Error: impossible to open 'data/points.txt'"<<std::endl;
      return false;
    }
    std::istream_iterator < Point > begin (in), end;
    T.insert (begin, end);
    CGAL::import_from_triangulation_3<LCC>(lcc,T);
    // Pb: the triangulation_3 is not the same on different machines ?
    // if ( !check_number_of_cells_3(lcc, 795, 4156, 6722, 3361, 1) )
    if ( !lcc.is_valid() )
      return false;
    lcc.clear();
  }
  
  return true;
}
