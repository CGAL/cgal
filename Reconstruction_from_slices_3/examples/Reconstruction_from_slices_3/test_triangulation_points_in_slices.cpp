// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : test_triangulation_points_in_slices.C
// package       : Reconstruction_from_slices
// author(s)     : Raphaelle Chaine
//
// ======================================================================

#include <iostream>

#ifdef GEOMVIEW
#include <CGAL/IO/Geomview_stream.h>
#endif //GEOMVIEW

#include <CGAL/Triangulation_from_slices_3.h>
#include <CGAL/Parser_SLC.h>

int main (int argc , char ** argv)
{
  const char * command = "Command :\n ./triangulation_points_in_slices slc_file [-geomview]";
  bool geomview    = false;

  // options of the command line
  if(argc == 1)
    {
      std::cout << command << std::endl;
      exit(1);
    }
  else if(argc > 3)
    {
      std::cout << command << std::endl;
      exit(2);
    }

  const char * slc_file = argv[1];
  for (int i=2 ; i<argc ; i++)
    {
      if (std::strcmp(argv[i],"-geomview")==0)
	geomview = true;
      else
	{
	  std::cout << command << std::endl;
	  exit(3);
	}
    }

  // Read data and triangulation construction
  typedef CGAL::Triangulation_from_slices_3< > Triangulation_from_slices_3;
  typedef CGAL::parser_SLC<Triangulation_from_slices_3> parser_SLC;

  Triangulation_from_slices_3 tpis;

  tpis.load<parser_SLC>(slc_file);

#ifdef GEOMVIEW
  if(geomview)
    {
      CGAL::Geomview_stream gv;
      gv.clear();
      gv << tpis;
      gv.look_recenter();
      std::getchar();
    }
#endif //GEOMVIEW

  return 0;
}


