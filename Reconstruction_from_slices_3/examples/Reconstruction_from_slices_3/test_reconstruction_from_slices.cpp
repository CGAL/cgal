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
// file          : test_reconstruction_from_slices.C
// package       : Nuage
// author(s)     : Raphaelle Chaine (Raphaelle.Chaine@sophia.inria.fr, raphaelle.chaine@liris.cnrs.fr)
//
// ======================================================================

#ifdef GEOMVIEW_DUMP //TO RM
#define GEOMVIEW     //TO RM
#include <CGAL/IO/Geomview_stream.h> //TO RM
CGAL::Geomview_stream gv;  //TO RM
#endif // GEOMVIEW_DUMP TO RM

#include <cstring> //strcat, strcpy
#include <cstdlib> //system
#include <iostream>

#ifdef GEOMVIEW //TO RM
#include <CGAL/IO/Geomview_stream.h> //TO RM
#endif //GEOMVIEW //TO RM

#include <CGAL/Triangulation_from_slices_3.h>
#include <CGAL/TFS_polyline_vertex_base_3.h>
#include <CGAL/TFS_cell_base_3.h>
#include <CGAL/Parser_CNT.h>
#include <CGAL/Reconstruction_from_slices_3.h>

int main (int argc , char ** argv)
{
#ifdef GEOMVIEW //TO RM
  const char * command = "Command :\n ./test_reconstruction_from_polylines cnt_file off_file_prefix \
[-solidlight or -solidheavy] [-heavy]  [-geomview]";
  bool geomview = false;
#else
  const char * command = "Command :\n ./test_reconstruction_from_polylines cnt_file off_file_prefix \
[-solidlight or -solidheavy] [-heavy]";
#endif //GEOMVIEW
  bool solidlight = false;
  bool solidheavy = false;
  bool heavy = false;

  // options of the command line
  if(argc == 1)
    {
      std::cout << command << std::endl;
      exit(1);
    }
  else if(argc > 6)
    {
      std::cout << command << std::endl;
      exit(2);
    }

  const char * cnt_file = argv[1];
  const char * off_file_prefix = argv[2];

  int l_pref=std::strlen(off_file_prefix);
  char *fname_off=new char[l_pref+10];
  std::strcat(std::strcpy(fname_off,off_file_prefix),".off");

  for (int i=3 ; i<argc ; i++)
    {
      if (std::strcmp(argv[i],"-solidlight")==0)
	solidlight = true;
      else if (std::strcmp(argv[i],"-solidheavy")==0)
	solidheavy = true;
      else if (std::strcmp(argv[i],"-heavy")==0)
	heavy = true;
#ifdef GEOMVIEW  //TO RM
      else if (std::strcmp(argv[i],"-geomview")==0)
	geomview = true;
#endif
      else
	{
	  std::cout << command << std::endl;
	  exit(3);
	}
    }

  // Read data and triangulation construction
  typedef CGAL::Exact_predicates_exact_constructions_kernel Gt;
  typedef CGAL::TFS_polyline_vertex_base_3< CGAL::Triangulation_vertex_base_3<Gt> > PlVb;
  typedef CGAL::TFS_cell_base_3< CGAL::Triangulation_cell_base_3<Gt> > PlCb;
  typedef CGAL::Triangulation_data_structure_3<PlVb,PlCb> Tds;
  typedef CGAL::Triangulation_from_slices_3< Gt,Tds> Slice_constrained_triangulation;

  Slice_constrained_triangulation tpis;

  typedef CGAL::parser_CNT<Slice_constrained_triangulation> CNT_parser;

  //   parser_CNT pp(cnt_file,&tpis);
  //   pp.parse();
  
  tpis.load<CNT_parser>(cnt_file);
  std::cout << "Triangulation over" << std::endl;

  if(!check_polylines_respect(tpis))
    {
      std::cout << "Not conform" << std::endl;
      check_and_conform_polylines(tpis);  
      std::cout << "Conforming over (guaranteed only for small angles)" << std::endl;
      CGAL_triangulation_assertion(check_polylines_respect(tpis));
    }
  {
    tag_inner_outer(tpis);
    std::cout << "In/out tagging over" << std::endl;
    if(solidlight)
      {
	non_solid_connections_removal(tpis,false); // Seems to be more efficient
	std::cout << "Non solid cluster removed." << std::endl;  
	std::cout << "Only edge-based local non solid connections remain." << std::endl;
      }
    else if(solidheavy)
      {
	non_solid_connections_removal(tpis,true);
	std::cout << "Non solid cluster removed." << std::endl; 
	std::cout << "Some edge-based local non solid connections are also removed (but some can remain)" << std::endl;
      }
    if(heavy)
      {
	stable_non_solid_connections_heavy_removal(tpis);
	std::cout << "Non solid cluster removed." << std::endl;	
        std::cout << "Edge-based local non solid connections removed." << std::endl; 
      }
    save_in_off_file(tpis,fname_off);
    std::cout << "Off file over" << std::endl;
    save_in_wrl_file(tpis, "hello.wrl",true);
  }
  delete fname_off;
  // else
  //   std::cout << "Non conformal polylines" << std::endl;
  
#ifdef GEOMVIEW
  if(geomview)
    {
      CGAL::Geomview_stream gvv;
      typedef Slice_constrained_triangulation::Finite_edges_iterator Finite_edges_iterator;
      typedef Slice_constrained_triangulation::Finite_cells_iterator Finite_cells_iterator;
      typedef Slice_constrained_triangulation::Segment Segment;
      typedef Slice_constrained_triangulation::Facet Facet;
      typedef Slice_constrained_triangulation::Triangle Triangle;
      typedef Slice_constrained_triangulation::Cell_handle Cell_handle;
      typedef Slice_constrained_triangulation::Cell_iterator Cell_iterator;
      gvv.clear();
      
      gvv.set_edge_color (CGAL::Color(0,0,0));
      gvv.set_face_color (CGAL::Color(0,0,0));
      gvv.set_bg_color (CGAL::Color(180,180,180));
      
      for (Cell_iterator it=tpis.cells_begin(); it != tpis.cells_end(); ++it)
	{
	  Cell_handle c(it);
	  if(c->is_external())
	    for(int i=0;i<4;i++)
	      if(c->neighbor(i)->is_internal()&&!c->neighbor(i)->is_non_solid())
		{
		  //gvv << Triangle(c->vertex((i+1)%4)->point(),c->vertex((i+2)%4)->point(),c->vertex((i+3)%4)->point());
		  gvv << Segment(c->vertex((i+1)%4)->point(),c->vertex((i+2)%4)->point());
		  gvv << Segment(c->vertex((i+2)%4)->point(),c->vertex((i+3)%4)->point());
		  gvv << Segment(c->vertex((i+3)%4)->point(),c->vertex((i+1)%4)->point());
		}
	}
      gvv.look_recenter();
      std::getchar();
    }
#endif //GEOMVIEW
  return 0;
}
