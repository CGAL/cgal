#include <CGAL/Nef_3/Single_wall_creator.h>
#include <CGAL/Nef_3/Single_wall_creator2.h>
#include <CGAL/Nef_3/Reflex_edge_searcher.h>
#include <CGAL/Nef_3/YVertical_wall_builder.h>
#include <CGAL/Nef_3/Ray_hit_generator2.h>
#include <CGAL/Nef_3/External_structure_builder.h>
#include <CGAL/Nef_3/SFace_separator.h>

#include <CGAL/Nef_3/Edge_sorter.h>
#include <CGAL/Nef_3/Edge_sorter2.h>

template<typename Nef_polyhedron>
void convex_decomposition_3(Nef_polyhedron& N) {

  typedef typename Nef_polyhedron::SNC_structure  SNC_structure;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename Nef_polyhedron::Vector_3           Vector_3;
  typedef typename Nef_polyhedron::Sphere_point   Sphere_point;

  typedef typename CGAL::Single_wall_creator<Nef_polyhedron>  Single_wall;
  typedef typename CGAL::Single_wall_creator2<Nef_polyhedron> Single_wall2;
  typedef typename CGAL::YVertical_wall_builder<Nef_polyhedron> YVertical_wall_builder;
  typedef typename YVertical_wall_builder::Vertical_redge_iterator Vertical_redge_iterator;
  typedef typename CGAL::Reflex_edge_searcher<Nef_polyhedron>  Reflex_edge_searcher;
  typedef typename Reflex_edge_searcher::Reflex_sedge_iterator Reflex_sedge_iterator;
  typedef typename Reflex_edge_searcher::Container            Container;
  typedef typename CGAL::Ray_hit_generator2<Nef_polyhedron> Ray_hit2;
  typedef typename CGAL::External_structure_builder<Nef_polyhedron> External_structure_builder;
  typedef typename CGAL::SFace_separator<Nef_polyhedron> SFace_separator;
  typedef typename CGAL::Edge_sorter<Nef_polyhedron, Container> Edge_sorter;
  typedef typename CGAL::Edge_sorter2<Nef_polyhedron, Container> Edge_sorter2;

  SFace_separator sf_sep;
  N.delegate(sf_sep,false, false);

  External_structure_builder esb;

  Reflex_edge_searcher res(Sphere_point(1,0,0));
  N.delegate(res,false,false);

  /*
  std::cerr << "number of reflex sedges" 
            << std::distance(res.positive_rsedges_begin(), 
                             res.positive_rsedges_end())
            << ","
            << std::distance(res.negative_rsedges_begin(), 
                             res.negative_rsedges_end())
            << std::endl;
  */

  Reflex_sedge_iterator rei;

  while(res.negative_rsedges_begin() != res.negative_rsedges_end()) {

  for(rei=res.negative_rsedges_begin(); rei!=res.negative_rsedges_end(); ++rei) {
      SHalfedge_handle se(*rei);
      Halfedge_handle split_edge;
    Ray_hit2 rh2a(Vector_3(-1,0,0),se->source()->source());
    N.delegate(rh2a);
    if(rh2a.split_edge(split_edge))
      res.handle_new_edge(split_edge);
    Ray_hit2 rh2b(Vector_3(-1,0,0),se->source()->twin()->source());
    N.delegate(rh2b);
    if(rh2b.split_edge(split_edge))
      res.handle_new_edge(split_edge);
  }
  
  /*
  int argc=0;
  char* argv[1];
  QApplication b(argc, argv);
  CGAL::Qt_widget_Nef_3<Nef_polyhedron>* wb = 
    new CGAL::Qt_widget_Nef_3<Nef_polyhedron>(N);
  b.setMainWidget(wb);
  wb->show();
  b.exec();
  */

  Edge_sorter es(res.get_negative_rsedges());
  N.delegate(es);

  //  CGAL_NEF_SETDTHREAD(227*229*47);
//  CGAL_NEF_SETDTHREAD(233*229*227);
  for(rei=res.negative_rsedges_begin(); rei!=res.negative_rsedges_end(); ++rei) {
    Halfedge_handle e = (*rei)->source();
//    std::cerr << "handle reflex edge " << e->source()->point() << "->" 
//    	      << e->twin()->source()->point() << std::endl;
    CGAL_assertion(res.is_reflex_sedge(*rei));    
    if(e->point().hx() > 0)
      e = e->twin();
    Single_wall W(e,Vector_3(-1,0,0));
    N.delegate(W);
  }    

  N.delegate(esb);

//  Reflex_edge_searcher res2;
  N.delegate(res, false, false);

  /*
  std::cerr << "number of reflex sedges" 
            << std::distance(res.positive_rsedges_begin(), 
                             res.positive_rsedges_end())
            << ","
            << std::distance(res.negative_rsedges_begin(), 
                             res.negative_rsedges_end())
            << std::endl;
  */
  }

  while(res.positive_rsedges_begin() != res.positive_rsedges_end()) {  
  for(rei=res.positive_rsedges_begin(); rei!=res.positive_rsedges_end(); ++rei) {
      SHalfedge_handle se(*rei);
      Halfedge_handle split_edge;
    Ray_hit2 rh2a(Vector_3(1,0,0),se->source()->source());
    N.delegate(rh2a);
    if(rh2a.split_edge(split_edge))
      res.handle_new_edge(split_edge);
    Ray_hit2 rh2b(Vector_3(1,0,0),se->source()->twin()->source());
    N.delegate(rh2b);
    if(rh2b.split_edge(split_edge))
      res.handle_new_edge(split_edge);
  }

  Edge_sorter2 es2(res.get_positive_rsedges());
  N.delegate(es2);


  //  CGAL_NEF_SETDTHREAD(229*233);

  rei=res.positive_rsedges_end();
  if(rei!=res.positive_rsedges_begin())
    do {
      --rei;
      Halfedge_handle e = (*rei)->source();
      //      std::cerr << "handle reflex edge " << e->source()->point() << "->" 
		//		<< e->twin()->source()->point() << std::endl;
      if(e->point().hx() < 0)
	e = e->twin();
      //      std::cerr << "handle reflex edge " << e->source()->point() << "->" 
      //		<< e->twin()->source()->point() << std::endl;
      Single_wall W(e,Vector_3(1,0,0));
      N.delegate(W);
    } while(rei!=res.positive_rsedges_begin());

  //  std::cerr << "here " << std::endl << N;

  N.delegate(esb);

  N.delegate(res, false, false);
  /*
  std::cerr << "number of reflex sedges" 
            << std::distance(res.positive_rsedges_begin(), 
                             res.positive_rsedges_end())
            << ","
            << std::distance(res.negative_rsedges_begin(), 
                             res.negative_rsedges_end())
            << std::endl;
  */

  } // TODO: get rid of the while loops

  YVertical_wall_builder Y;
  N.delegate(Y,false,false);

  typename YVertical_wall_builder::Vertical_redge_iterator vri=Y.pos_begin();
  for(; vri != Y.pos_end(); ++vri) {
    //    std::cerr << "pos: " << (*vri)->source()->point()
    //	      << "->" << (*vri)->twin()->source()->point() << std::endl;
    Single_wall2 W((*vri),Sphere_point(0,1,0));
    N.delegate(W);
  }
  for(vri = Y.neg_begin(); vri != Y.neg_end(); ++vri) {
    //    std::cerr << "neg: " << (*vri)->source()->point()
    //	      << "->" << (*vri)->twin()->source()->point() << std::endl;
    Single_wall2 W((*vri),Sphere_point(0,-1,0));
    N.delegate(W);
  }      

  N.delegate(esb);
}
