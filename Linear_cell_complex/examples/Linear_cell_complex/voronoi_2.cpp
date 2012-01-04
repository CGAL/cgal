// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <iostream>
#include <fstream>

/* If you want to use a viewer, you can use one of the following file
 * depending if you use vtk or qglviewer. */
#ifdef CGAL_LCC_USE_QT
#include "linear_cell_complex_3_viewer_qt.h"
#else 
#ifdef CGAL_LCC_USE_VTK
#include "linear_cell_complex_3_viewer_vtk.h"
#endif
#endif

/* // If you want to use exact constructions.
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
typedef CGAL::Linear_cell_complex<2,2,
  CGAL::Linear_cell_complex_traits<2, CGAL::Exact_predicates_exact_constructions_kernel> > LCC_2;
*/

typedef CGAL::Linear_cell_complex<2> LCC_2;
typedef LCC_2::Dart_handle           Dart_handle;
typedef LCC_2::Point                 Point;

typedef CGAL::Delaunay_triangulation_2<LCC_2::Traits> Triangulation;

// Function used to display the voronoi diagram.
void display_voronoi(LCC_2& alcc, Dart_handle adart)
{
  // We remove the infinite face plus all the faces adjacent to it.
  // Indeed, we cannot view these faces since they do not have
  // a "correct geometry". 
  std::stack<Dart_handle> toremove;
  int mark_toremove=alcc.get_new_mark();

  // adart belongs to the infinite face.
  toremove.push(adart);
  CGAL::mark_cell<LCC_2,2>(alcc, adart, mark_toremove);

  // Now we get all the faces adjacent to the infinite face.
  for (LCC_2::Dart_of_cell_range<2>::iterator
         it=alcc.darts_of_cell<2>(adart).begin(),
         itend=alcc.darts_of_cell<2>(adart).end(); it!=itend; ++it)
  {
    if ( !alcc.is_marked(it->beta(2), mark_toremove) )
    {
      CGAL::mark_cell<LCC_2,2>(alcc, it->beta(2), mark_toremove);
      toremove.push(it->beta(2));
    }    
  }
  
  while( !toremove.empty() )
  {
    CGAL::remove_cell<LCC_2, 2>(alcc, toremove.top());
    toremove.pop();
  }

  CGAL_assertion(alcc.is_without_boundary(1));
    
  std::cout<<"Voronoi subdvision, only finite faces:"<<std::endl<<"  ";
  alcc.display_characteristics(std::cout) << ", valid=" 
                                          << alcc.is_valid()
                                          << std::endl;
  
#ifdef CGAL_LCC_USE_VIEWER
  display_lcc(alcc);
#endif // CGAL_LCC_USE_VIEWER
}

template<typename LCC, typename TR>
void transform_dart_to_their_dual(LCC& alcc, LCC& adual,
                                  std::map<typename TR::Face_handle,
                                           typename LCC::Dart_handle>& assoc)
{
  typename LCC::Dart_range::iterator it1=alcc.darts().begin();
  typename LCC::Dart_range::iterator it2=adual.darts().begin();

  std::map<typename LCC::Dart_handle, typename LCC::Dart_handle> dual;
  
  for ( ; it1!=alcc.darts().end(); ++it1, ++it2 )
  {
    dual[it1]=it2;
  }

  for ( typename std::map<typename TR::Face_handle, typename LCC::Dart_handle>
          ::iterator it=assoc.begin(), itend=assoc.end(); it!=itend; ++it)
  {
    assoc[it->first]=dual[it->second];
  } 
}

template<typename LCC, typename TR>
void set_geometry_of_dual(LCC& alcc, TR& tr,
                          std::map<typename TR::Face_handle,
                                   typename LCC::Dart_handle>& assoc)
{
  for ( typename std::map<typename TR::Face_handle, typename LCC::Dart_handle>
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
    filename=std::string("data/points_2");
    std::cout<<"No filename given: use data/points_2 by default."<<std::endl;
  }
  else
    filename=std::string(argv[1]);
  
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
           LCC_2::Dart_handle > face_to_dart;
  
  Dart_handle dh=CGAL::import_from_triangulation_2<LCC_2, Triangulation>
    (lcc, T, &face_to_dart);
  CGAL_assertion(lcc.is_without_boundary());

  std::cout<<"Delaunay triangulation :"<<std::endl<<"  ";
  lcc.display_characteristics(std::cout) << ", valid=" 
                                         << lcc.is_valid() << std::endl;

  // 3) Compute the dual lcc.
  LCC_2 dual_lcc;
  Dart_handle ddh=lcc.dual(dual_lcc, dh);
  // Here, dual_lcc is the 2D Voronoi diagram.
  CGAL_assertion(dual_lcc.is_without_boundary());

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

  display_voronoi(dual_lcc, ddh);

  return EXIT_SUCCESS;
}

