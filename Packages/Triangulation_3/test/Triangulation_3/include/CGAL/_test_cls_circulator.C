// ============================================================================
//
// Copyright (c) 1998 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// source        : 
// file          : include/CGAL/_test_cls_circulator.C
// revision      : 
// revision_date : 
// author(s)     : Francois Rebufat (Francois.Rebufat@sophia.inria.fr)
//
// coordinator   : INRIA Sophia-Antipolis
// ============================================================================



template < class Triangulation >
int
_test_circulator( const Triangulation &T )
{
  typedef typename Triangulation::Edge_iterator    Edge_iterator;
  typedef typename Triangulation::Vertex_iterator  Vertex_iterator;
  typedef typename Triangulation::Cell_circulator  Cell_circulator;
  typedef typename Triangulation::Cell             Cell;
  typedef typename Triangulation::Vertex           Vertex;
  typedef typename Triangulation::Vertex_handle    Vertex_handle;
  typedef typename Triangulation::Cell_handle      Cell_handle;
  //  typedef typename Triangulation::Facet_circulator Facet_circulator;


  int n = 0;
  Cell_circulator cc, cc0;
  Edge_iterator eit;
  // testing incident_cells(edge *);
  for (eit=T.all_edges_begin(); eit!=T.edges_end();  eit++)
    {
     cc0=cc=T.incident_cells(*eit);
      do {
	cc++; n++;
      } while (cc != cc0);
    }
  // testing incident_cells(edge *,cell *);
  for (eit=T.all_edges_begin(); eit!=T.edges_end();  eit++)
    {
     cc0=cc=T.incident_cells(*eit,(*eit).first);
      do {
	cc++; n++;
      } while (cc != cc0);
    }
  for (eit=T.finite_edges_begin(); eit!=T.edges_end();  eit++)
    {
     cc0=cc=T.incident_cells(*eit);
      do {
	cc++; n++;
      } while (cc != cc0);
    }
  for (eit=T.finite_edges_begin(); eit!=T.edges_end(); eit++)
    {
     cc0=cc=T.incident_cells(*eit,(*eit).first);
      do {
	cc++; n++;
      } while (cc != cc0);
    }
  std::set<Cell*, less<Cell*> > cells ;
  std::set<Vertex*, less<Vertex*> > vertices ;

  Vertex_iterator vit;
  for (vit=T.all_vertices_begin(); vit!=T.vertices_end() ; vit++) {
    T.incident_cells(vit,cells);
    T.incident_cells(vit,cells,vit->cell());
    T.incident_vertices(vit, vertices);
    T.incident_vertices(vit, vertices,vit->cell());
    
  }

//    Facet_circulator fc, fc0;
//    for (eit=T.all_edges_begin(); eit!=T.edges_end(); eit++)
//     {
//      fc0=fc=T.incident_facets(*eit);
//       do {
// 	fc++; n++;
//       } while (fc != fc0);
//     }
//    for (eit=T.all_edges_begin(); eit!=T.edges_end(); eit++)
//     {
//      fc0=fc=T.incident_facets(*eit,eit-->facet());
//       do {
// 	fc++; n++;
//       } while (fc != fc0);
//     }
//    int fi;
//    for (eit=T.all_edges_begin(); eit!=T.edges_end(); eit++)
//     {
//      for (fi=0; fi!=4 ; fi++) 
//        {
// 	if (t.dimension()==2) {fi=3;}
//         fc0=fc=T.incident_facets(*eit,eit->cell(),fi);
//         do {
//    	  fc++; n++;
//         } while (fc != fc0);
//        }
//     }

//    for (eit=T.finite_edges_begin(); eit!=T.edges_end(); eit++)
//     {
//      fc0=fc=T.incident_facets(*eit);
//       do {
// 	fc++; n++;
//       } while (fc != fc0);
//     }
//    for (eit=T.finite_edges_begin(); eit!=T.edges_end(); eit++)
//     {
//      fc0=fc=T.incident_facets(*eit,eit-->facet());
//       do {
// 	fc++; n++;
//       } while (fc != fc0);
//     }
//    int fi;
//    for (eit=T.finite_edges_begin(); eit!=T.edges_end(); eit++)
//     {
//      for (fi=0; fi!=4 ; fi++) 
//        {
// 	if (t.dimension()==2) {fi=3;}
//         fc0=fc=T.incident_facets(*eit,eit->cell(),fi);
//         do {
//    	  fc++; n++;
//         } while (fc != fc0);
//        }
//     }
  return n;
}
