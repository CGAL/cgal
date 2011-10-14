// Copyright (c) 2010 CNRS, LIRIS, http://liris.cnrs.fr/, All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
#include <CGAL/Timer.h>
#include <iostream>
#include <fstream>

#define NUMBER_OF_LOOP 10

typedef CGAL::Combinatorial_map_with_points<3> CMap_3;
typedef CMap_3::Dart_handle          Dart_handle;

int main(int narg, char** argv)
{
  if ( narg==1 ||
       (narg>1 && (!strcmp(argv[1],"-h") || !strcmp(argv[1],"-?"))) )
    {
      std::cout<<"Usage : a.out [-h -?] filename1 ... filenamek"<<std::endl;
      std::cout<<"  Load all the filename into a same map "
	"each file must be an off file) and run iterators tests."<<std::endl;
      exit(EXIT_SUCCESS);
    }

  CMap_3 m;
  
  for (int i=1; i<narg; ++i)
    {
      std::ifstream is(argv[i]);
      CGAL::import_from_polyhedron_flux<CMap_3>(m,is);
      is.close();
    }

  CMap_3::size_type res=0;
  CGAL::Timer  t; 

  std::cout << "*****************************************" << std::endl;
  std::cout << "****** Iterator based on Container ******" << std::endl
	    << std::endl;

  std::cout << "****** Iterators with range ******" << std::endl;
  res=0;
  t.reset();
  t.start();
  //  for (unsigned int i=0; i<NUMBER_OF_LOOP; ++i)
  for (CMap_3::Dart_range::iterator it(m.darts().begin()), itend(m.darts().end()); it!=itend; ++it)
    {
      for (CMap_3::Vertex_attribute_range::iterator it2(m.vertex_attributes().begin()); 
	   it2!=m.vertex_attributes().end(); ++it2)
	{
	  ++res;
	}
    }
  t.stop();
  std::cout<<"Result: "<<res<<" in "<<t.time()<<" seconds."<<std::endl;

  {
    CMap_3::Dart_of_orbit_range<2>::iterator 
      it1=m.darts_of_orbit<2>(m.first_dart()).begin();
    CMap_3::Dart_of_orbit_range<2>::const_iterator it2=it1;
  }
  {
    CMap_3::Dart_of_involution_range<2>::iterator 
      it1(m,m.first_dart());
    CMap_3::Dart_of_involution_range<2>::const_iterator it2=it1;
  }
  {
    CMap_3::Dart_of_cell_range<2>::iterator 
      it1=m.darts_of_cell<2>(m.first_dart()).begin(); 
    CMap_3::Dart_of_cell_range<2>::const_iterator it2=it1;
  }
  {
    CMap_3::One_dart_per_incident_cell_range<2,0>::iterator 
      it1(m,m.first_dart());
    CMap_3::One_dart_per_incident_cell_range<2,0>::const_iterator it2=it1;
  }

  std::cout << "****** Iterators with range, with end declare before the loop ******" << std::endl;
  res=0;
  t.reset();
  t.start();
  {
    CMap_3::Dart_range::iterator itend(m.darts().end());
    CMap_3::Vertex_attribute_range::iterator it2end(m.vertex_attributes().end());
    //  for (unsigned int i=0; i<NUMBER_OF_LOOP; ++i)
    for (CMap_3::Dart_range::iterator it(m.darts().begin()); it!=itend; ++it)
      for (CMap_3::Vertex_attribute_range::iterator it2(m.vertex_attributes().begin()); 
	   it2!=it2end; ++it2)
	{
	  ++res;
	}
  }
  t.stop();
  std::cout<<"Result: "<<res<<" in "<<t.time()<<" seconds."<<std::endl;

  std::cout << std::endl << "****************************" << std::endl;
  std::cout << "****** Other iterator ******" << std::endl;

  std::cout << "****** Iterators with cont ******" << std::endl;
  res=0;
  t.reset();
  t.start();
  for (unsigned int i=0; i<NUMBER_OF_LOOP; ++i)
    for (CMap_3::Dart_range::iterator it(m.darts().begin()); 
	 it!=m.darts().end(); ++it)
      {
	for (CGAL::CMap_one_dart_per_incident_cell_iterator<CMap_3,2,0> it2(m,it); 
	     it2.cont(); ++it2)
	  ++res;
      }
  t.stop();
  std::cout<<"Result: "<<res<<" in "<<t.time()<<" seconds."<<std::endl;

  std::cout << "****** Iterators with range ******" << std::endl;
  res=0;
  t.reset();
  t.start();
  for (unsigned int i=0; i<NUMBER_OF_LOOP; ++i)
    for (CMap_3::Dart_range::iterator it(m.darts().begin()); 
	 it!=m.darts().end(); ++it)
      {
	for (CMap_3::One_dart_per_incident_cell_range<2,0>::iterator 
	       it2=m.one_dart_per_incident_cell<2,0>(it).begin();
	     it2!=m.one_dart_per_incident_cell<2,0>(it).end(); ++it2)
	  ++res;
      }
  t.stop();
  std::cout<<"Result: "<<res<<") in "<<t.time()<<" seconds."<<std::endl;

  std::cout << "****** Iterators with range and end, with end declare before the loop ******" << std::endl;
  res=0;
  t.reset();
  t.start();
  {
    CMap_3::Dart_range::iterator itend(m.darts().end());
    CMap_3::One_dart_per_incident_cell_range<2,0>::iterator 
      it2end(m.one_dart_per_incident_cell<2,0>(NULL).end());
    for (unsigned int i=0; i<NUMBER_OF_LOOP; ++i)
      {
	for (CMap_3::Dart_range::iterator it(m.darts().begin()); it!=itend; ++it)
	  {
	    for (CMap_3::One_dart_per_incident_cell_range<2,0>::iterator
		   it2(m.one_dart_per_incident_cell<2,0>(it).begin());
		 it2!=it2end; ++it2)
	      ++res;
	  }
      }
  }
  t.stop();
  std::cout<<"Result: "<<res<<" in "<<t.time()<<" seconds."<<std::endl;

  std::cout << "****** Iterators with range 'Sylvain' version ******" << std::endl;
  res=0;
  t.reset();
  t.start();
  for (unsigned int i=0; i<NUMBER_OF_LOOP; ++i)
    for (CMap_3::Dart_range::iterator it(m.darts().begin()), 
	   itend(m.darts().end()); it!=itend; ++it)
      {
	for (CMap_3::One_dart_per_incident_cell_range<2,0>::iterator 
	       it2=m.one_dart_per_incident_cell<2,0>(it).begin(),
	       it2end(m.one_dart_per_incident_cell<2,0>(it).end()); 
	     it2!=it2end; ++it2)
	  ++res;
      }
  t.stop();
  std::cout<<"Result: "<<res<<" in "<<t.time()<<" seconds."<<std::endl;

  return 1;
}
