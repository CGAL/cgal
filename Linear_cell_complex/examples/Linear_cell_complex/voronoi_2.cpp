#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_2_to_lcc.h>

#include <iostream>
#include <fstream>
#include <cassert>


// This example works both with cmap and gmap as combinatorial data structure.
//typedef CGAL::Linear_cell_complex_for_combinatorial_map<2> LCC_2;
typedef CGAL::Linear_cell_complex_for_generalized_map<2> LCC_2;
typedef LCC_2::Dart_descriptor Dart_descriptor;
typedef LCC_2::Point           Point;

typedef CGAL::Delaunay_triangulation_2<LCC_2::Traits> Triangulation;

// Function used to display the voronoi diagram.
void display_voronoi(LCC_2& alcc, Dart_descriptor adart)
{
  // We remove the infinite face plus all the faces adjacent to it.
  // Indeed, we cannot view these faces since they do not have
  // a "correct geometry".
  std::stack<Dart_descriptor> toremove;
  LCC_2::size_type mark_toremove=alcc.get_new_mark();

  // adart belongs to the infinite face.
  toremove.push(adart);
  CGAL::mark_cell<LCC_2,2>(alcc, adart, mark_toremove);

  // Now we get all the faces adjacent to the infinite face.
  for (LCC_2::Dart_of_cell_range<2>::iterator
         it=alcc.darts_of_cell<2>(adart).begin(),
         itend=alcc.darts_of_cell<2>(adart).end(); it!=itend; ++it)
  {
    if ( !alcc.is_marked(alcc.opposite<2>(it), mark_toremove) )
    {
      CGAL::mark_cell<LCC_2,2>(alcc, alcc.opposite<2>(it), mark_toremove);
      toremove.push(alcc.opposite<2>(it));
    }
  }

  while( !toremove.empty() )
  {
    alcc.remove_cell<2>(toremove.top());
    toremove.pop();
  }

  assert(alcc.is_without_boundary(1));

  std::cout<<"Voronoi subdvision, only finite faces:"<<std::endl<<"  ";
  alcc.display_characteristics(std::cout) << ", valid="
                                          << alcc.is_valid()
                                          << std::endl;
}

template<typename LCC, typename TR>
void transform_dart_to_their_dual(LCC& alcc, LCC& adual,
                                  std::map<typename TR::Face_handle,
                                           typename LCC::Dart_descriptor>& assoc)
{
  typename LCC::Dart_range::iterator it1=alcc.darts().begin();
  typename LCC::Dart_range::iterator it2=adual.darts().begin();

  std::map<typename LCC::Dart_descriptor, typename LCC::Dart_descriptor> dual;

  for ( ; it1!=alcc.darts().end(); ++it1, ++it2 )
  {
    dual[it1]=it2;
  }

  for ( typename std::map<typename TR::Face_handle, typename LCC::Dart_descriptor>
          ::iterator it=assoc.begin(), itend=assoc.end(); it!=itend; ++it)
  {
    assoc[it->first]=dual[it->second];
  }
}

template<typename LCC, typename TR>
void set_geometry_of_dual(LCC& alcc, TR& tr,
                          std::map<typename TR::Face_handle,
                                   typename LCC::Dart_descriptor>& assoc)
{
  for ( typename std::map<typename TR::Face_handle, typename LCC::Dart_descriptor>
          ::iterator it=assoc.begin(), itend=assoc.end(); it!=itend; ++it)
  {
    if ( !tr.is_infinite(it->first) )
      alcc.set_vertex_attribute
        (it->second,alcc.create_vertex_attribute(tr.circumcenter(it->first)));
    else
      alcc.set_vertex_attribute(it->second,alcc.create_vertex_attribute());
  }
}

int main(int narg, char** argv)
{
  if (narg>1 && (!strcmp(argv[1],"-h") || !strcmp(argv[1],"-?")) )
  {
    std::cout<<"Usage : voronoi_2 filename"<<std::endl
             <<"   filename being a fine containing 2D points used to "
             <<" compute the Delaunay_triangulation_2."<<std::endl;
    return EXIT_FAILURE;
  }

  std::string filename;
  if ( narg==1 )
  {
    filename="data/points_2";
    std::cout<<"No filename given: use "<<filename<<" by default."<<std::endl;
  }
  else { filename=std::string(argv[1]); }

  // 1) Compute the Delaunay_triangulation_2.
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

  // 2) Convert the triangulation into a 2D lcc.
  LCC_2 lcc;
  std::map<Triangulation::Face_handle,
           LCC_2::Dart_descriptor > face_to_dart;

  Dart_descriptor d=CGAL::import_from_triangulation_2<LCC_2, Triangulation>
    (lcc, T, &face_to_dart);
  assert(lcc.is_without_boundary());

  std::cout<<"Delaunay triangulation :"<<std::endl<<"  ";
  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid() << std::endl;

  // 3) Compute the dual lcc.
  LCC_2 dual_lcc;
  Dart_descriptor dd=lcc.dual(dual_lcc, d);
  // Here, dual_lcc is the 2D Voronoi diagram.
  assert(dual_lcc.is_without_boundary());

  // 4) We update the geometry of dual_lcc by using the std::map
  //    face_to_dart.
  transform_dart_to_their_dual<LCC_2,Triangulation>
    (lcc, dual_lcc, face_to_dart);
  set_geometry_of_dual<LCC_2,Triangulation>(dual_lcc, T, face_to_dart);

  // 5) Display the dual_lcc characteristics.
  std::cout<<"Voronoi subdvision :"<<std::endl<<"  ";
  dual_lcc.display_characteristics(std::cout) << ", valid="
                                              << dual_lcc.is_valid()
                                              << std::endl;

  display_voronoi(dual_lcc, dd);

  return EXIT_SUCCESS;
}
