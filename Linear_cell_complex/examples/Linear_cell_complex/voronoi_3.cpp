#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Triangulation_3_to_lcc.h>

#include <iostream>
#include <fstream>
#include <cassert>

/* // If you want to use exact constructions.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Linear_cell_complex<3,3,
  CGAL::Linear_cell_complex_traits<3, CGAL::Exact_predicates_exact_constructions_kernel> > LCC_3;
*/

typedef CGAL::Linear_cell_complex_for_combinatorial_map<3> LCC_3;
typedef LCC_3::Dart_descriptor Dart_descriptor;
typedef LCC_3::Point           Point;

typedef CGAL::Delaunay_triangulation_3<LCC_3::Traits> Triangulation;

// Function used to display the voronoi diagram.
void display_voronoi(LCC_3& alcc, Dart_descriptor adart)
{
  // We remove the infinite volume plus all the volumes adjacent to it.
  // Indeed, we cannot view these volumes since they do not have
  // a "correct geometry".
  std::stack<Dart_descriptor> toremove;
  LCC_3::size_type mark_toremove=alcc.get_new_mark();

  // adart belongs to the infinite volume.
  toremove.push(adart);
  CGAL::mark_cell<LCC_3,3>(alcc, adart, mark_toremove);

  // Now we get all the volumes adjacent to the infinite volume.
  for (LCC_3::Dart_of_cell_range<3>::iterator
         it=alcc.darts_of_cell<3>(adart).begin(),
         itend=alcc.darts_of_cell<3>(adart).end(); it!=itend; ++it)
  {
    if ( !alcc.is_marked(alcc.beta(it,3), mark_toremove) )
    {
      CGAL::mark_cell<LCC_3,3>(alcc, alcc.beta(it,3), mark_toremove);
      toremove.push(alcc.beta(it,3));
    }
  }

  while( !toremove.empty() )
  {
    alcc.remove_cell<3>(toremove.top());
    toremove.pop();
  }

  assert(alcc.is_without_boundary(1) && alcc.is_without_boundary(2));

  std::cout<<"Voronoi subdvision, only finite volumes:"<<std::endl<<"  ";
  alcc.display_characteristics(std::cout) << ", valid="
                                          << alcc.is_valid()
                                          << std::endl;
}

template<typename LCC, typename TR>
void transform_dart_to_their_dual(LCC& alcc, LCC& adual,
                                  std::map<typename TR::Cell_handle,
                                           typename LCC::Dart_descriptor>& assoc)
{
  typename LCC::Dart_range::iterator it1=alcc.darts().begin();
  typename LCC::Dart_range::iterator it2=adual.darts().begin();

  std::map<typename LCC::Dart_descriptor, typename LCC::Dart_descriptor> dual;

  for ( ; it1!=alcc.darts().end(); ++it1, ++it2 )
  {
    dual[it1]=it2;
  }

  for ( typename std::map<typename TR::Cell_handle, typename LCC::Dart_descriptor>
          ::iterator it=assoc.begin(), itend=assoc.end(); it!=itend; ++it)
  {
    assoc[it->first]=dual[it->second];
  }
}

template<typename LCC, typename TR>
void set_geometry_of_dual(LCC& alcc, TR& tr,
                          std::map<typename TR::Cell_handle,
                                   typename LCC::Dart_descriptor>& assoc)
{
  for ( typename std::map<typename TR::Cell_handle, typename LCC::Dart_descriptor>
          ::iterator it=assoc.begin(), itend=assoc.end(); it!=itend; ++it)
  {
    if ( !tr.is_infinite(it->first) )
      alcc.set_vertex_attribute
        (it->second,alcc.create_vertex_attribute(tr.dual(it->first)));
    else
      alcc.set_vertex_attribute(it->second,alcc.create_vertex_attribute());
  }
}


int main(int narg, char** argv)
{
  if (narg>1 && (!strcmp(argv[1],"-h") || !strcmp(argv[1],"-?")) )
  {
    std::cout<<"Usage : voronoi_3 filename"<<std::endl
             <<"   filename being a fine containing 3D points used to "
             <<" compute the Delaunay_triangulation_3."<<std::endl;
    return EXIT_FAILURE;
  }

  std::string filename;
  if ( narg==1 )
  {
    filename=std::string("data/points_3");
    std::cout<<"No filename given: use data/points_3 by default."<<std::endl;
  }
  else
    filename=std::string(argv[1]);

  // 1) Compute the Delaunay_triangulation_3.
  Triangulation T;

  std::ifstream iFile(filename.c_str());
  if (!iFile)
  {
    std::cout << "Problem reading file " << filename << std::endl;
    return EXIT_FAILURE;
  }

  std::istream_iterator<Point> begin(iFile), end;
  T.insert(begin, end);
  assert(T.is_valid(false));

  // 2) Convert the triangulation into a 3D lcc.
  LCC_3 lcc;
  std::map<Triangulation::Cell_handle,
           LCC_3::Dart_descriptor > vol_to_dart;

  Dart_descriptor d=CGAL::import_from_triangulation_3<LCC_3, Triangulation>
    (lcc, T, &vol_to_dart);

  std::cout<<"Delaunay triangulation :"<<std::endl<<"  ";
  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid() << std::endl;

  // 3) Compute the dual lcc.
  LCC_3 dual_lcc;
  Dart_descriptor dd=lcc.dual(dual_lcc, d);
  // Here, dual_lcc is the 3D Voronoi diagram.
  assert(dual_lcc.is_without_boundary());

  // 4) We update the geometry of dual_lcc by using the std::map
  //    face_to_dart.
  transform_dart_to_their_dual<LCC_3,Triangulation>
    (lcc, dual_lcc, vol_to_dart);
  set_geometry_of_dual<LCC_3,Triangulation>(dual_lcc, T, vol_to_dart);

  // 5) Display the dual_lcc characteristics.
  std::cout<<"Voronoi subdvision :"<<std::endl<<"  ";
  dual_lcc.display_characteristics(std::cout) << ", valid="
                                              << dual_lcc.is_valid()
                                              << std::endl;
  display_voronoi(dual_lcc, dd);

  return EXIT_SUCCESS;
}
