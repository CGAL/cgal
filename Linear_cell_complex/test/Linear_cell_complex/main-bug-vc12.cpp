#include <CGAL/Linear_cell_complex.h>
#include <iostream>

template<typename LCC>
bool check_number_of_cells_2(LCC& lcc, unsigned int nbv, unsigned int nbe,
                             unsigned int nbf, unsigned int nbcc)
{
  if ( !lcc.is_valid() )
    {
      std::cout<<"ERROR: the lcc is not valid."<<std::endl;
      assert(false);
      return false;
    }

  std::vector<unsigned int> nbc;
  nbc=lcc.count_all_cells();

  if (nbv!=nbc[0] || nbe!=nbc[1] || nbf!=nbc[2] || nbcc!=nbc[3])
    {
      std::cout<<"ERROR: the number of cells is not correct. We must have "
               <<" ("<<nbv<<", "<<nbe<<", "<<nbf<<", "<<nbcc<<") and we have"
               <<" ("<<nbc[0]<<", "<<nbc[1]<<", "<<nbc[2]<<", "<<nbc[3]<<")."
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

  std::cout<<"OK"<<std::endl;

  return true;
}


int main()
{
  typedef CGAL::Linear_cell_complex<2> LCC;
  typedef LCC::Dart_handle Dart_handle;
  typedef LCC::Point Point;
  
  LCC lcc;
  
  std::cout<<"Make 2 triangles..."<<std::flush;
  Dart_handle dh1=lcc.make_triangle(Point(5,5),Point(7,5),Point(6,6));
  Dart_handle dh2=lcc.make_triangle(Point(5,4),Point(7,4),Point(6,3));
  if ( !check_number_of_cells_2(lcc, 6, 6, 2, 2) )
    return EXIT_FAILURE;

  std::cout<<"2-sew the 2 triangles..."<<std::flush;
  lcc.template sew<2>(dh1,dh2);
  if ( !check_number_of_cells_2(lcc, 4, 5, 2, 1) )
    return EXIT_FAILURE;

  std::cout<<"2-unsew the 2 triangles..."<<std::flush;
  lcc.template unsew<2>(dh1);
  if ( !check_number_of_cells_2(lcc, 6, 6, 2, 2) )
    return EXIT_FAILURE;
  
  return EXIT_SUCCESS;
}
