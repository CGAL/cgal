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
// file          : test_triangulation_points_in_polylines.C
// package       : Nuage
// author(s)     : Raphaelle Chaine  (Raphaelle.Chaine@sophia.inria.fr, raphaelle.chaine@liris.cnrs.fr)
//
// ======================================================================

#include <iostream>

#ifdef GEOMVIEW
#include <CGAL/IO/Geomview_stream.h>
#endif //GEOMVIEW

#include <CGAL/Triangulation_from_slices_3.h>
#include <CGAL/TFS_polyline_vertex_base_3.h>
#include <CGAL/Parser_CNT.h>

int main (int argc , char ** argv)
{
  const char * command = "Command :\n ./test_triangulation_points_in_polylines cnt_file [-geomview]";
  bool geomview = false;

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
  const char * cnt_file = argv[1];
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
  typedef CGAL::Exact_predicates_exact_constructions_kernel Gt;
  typedef CGAL::TFS_polyline_vertex_base_3<CGAL::Triangulation_vertex_base_3<Gt> > Vertex;
  typedef CGAL::Triangulation_cell_base_3<Gt> Cell;
  typedef CGAL::Triangulation_data_structure_3<Vertex,Cell> Tds;
  typedef CGAL::Triangulation_from_slices_3<Gt,Tds> Triangulation_from_slices_3;

  Triangulation_from_slices_3 tpis;

  typedef CGAL::parser_CNT<Triangulation_from_slices_3> parser_CNT;

  //   parser_CNT pp(cnt_file,&tpis);
  //   pp.parse();

  tpis.load<parser_CNT>(cnt_file);

#ifdef GEOMVIEW
  if(geomview)
    {
      typedef Triangulation_from_slices_3::Finite_edges_iterator Finite_edges_iterator;
      typedef Triangulation_from_slices_3::Segment Segment;

      CGAL::Geomview_stream gv;
      gv.clear();
      gv.set_edge_color (CGAL::Color(0,0,0));
      //gv.set_face_color (CGAL::Color(255,255,255));
      gv.set_bg_color (CGAL::Color(180,180,180));

      for (Finite_edges_iterator it=tpis.finite_edges_begin(); it != tpis.finite_edges_end(); ++it)
	if (((*it).first->vertex((*it).second)->get_next()==(*it).first->vertex((*it).third))
	    ||((*it).first->vertex((*it).third)->get_next()==(*it).first->vertex((*it).second)))
	  {
	    gv << Segment((*it).first->vertex((*it).second)->point(),(*it).first->vertex((*it).third)->point());
	  }
      gv.look_recenter();
      std::getchar();
    }
#endif //GEOMVIEW
  return 0;
}


