#include <CGAL/Nef_3/Single_wall_creator.h>
#include <CGAL/Nef_3/Single_wall_creator2.h>
#include <CGAL/Nef_3/Single_wall_creator3.h>
#include <CGAL/Nef_3/Reflex_edge_searcher.h>
#include <CGAL/Nef_3/YVertical_wall_builder.h>
#include <CGAL/Nef_3/YVertical_wall_builder2.h>
#include <CGAL/Nef_3/Ray_hit_generator2.h>
#include <CGAL/Nef_3/External_structure_builder.h>
#include <CGAL/Nef_3/SFace_separator.h>

#include <CGAL/Nef_3/Edge_sorter.h>
#include <CGAL/Nef_3/Edge_sorter2.h>

#include <CGAL/Nef_3/is_reflex_sedge.h>

template<typename Nef_polyhedron>
void convex_decomposition_3(Nef_polyhedron& N) {

  typedef typename Nef_polyhedron::SNC_structure  SNC_structure;
  typedef typename SNC_structure::Vertex_handle Vertex_handle;
  typedef typename SNC_structure::Halfedge_handle Halfedge_handle;
  typedef typename SNC_structure::SHalfedge_handle SHalfedge_handle;
  typedef typename SNC_structure::SHalfedge_iterator SHalfedge_iterator;
  typedef typename Nef_polyhedron::Point_3            Point_3;
  typedef typename Nef_polyhedron::Vector_3           Vector_3;
  typedef typename Nef_polyhedron::Sphere_point   Sphere_point;
  typedef typename Nef_polyhedron::FT FT;

  typedef typename CGAL::Single_wall_creator<Nef_polyhedron>  Single_wall;
  typedef typename CGAL::Single_wall_creator2<Nef_polyhedron> Single_wall2;
  typedef typename CGAL::Single_wall_creator3<Nef_polyhedron> Single_wall3;
  typedef typename CGAL::YVertical_wall_builder<Nef_polyhedron> YVertical_wall_builder;
  typedef typename CGAL::YVertical_wall_builder2<Nef_polyhedron> YVertical_wall_builder2;
  typedef typename CGAL::Reflex_vertex_searcher<Nef_polyhedron>  Reflex_vertex_searcher;
  typedef typename Reflex_vertex_searcher::Reflex_vertex_map Reflex_vertex_map;
  typedef typename CGAL::Reflex_edge_searcher<Nef_polyhedron, Reflex_vertex_map>  Reflex_edge_searcher;
  typedef typename Reflex_edge_searcher::Reflex_sedge_iterator Reflex_sedge_iterator;
  typedef typename Reflex_edge_searcher::Container            Container;
  typedef typename CGAL::Ray_hit_generator2<Nef_polyhedron> Ray_hit2;
  typedef typename CGAL::External_structure_builder<Nef_polyhedron> External_structure_builder;
  typedef typename CGAL::SFace_separator<Nef_polyhedron> SFace_separator;
  typedef typename CGAL::Edge_sorter<Nef_polyhedron, std::less<Point_3>, 
                                     std::less<FT>, Container> Edge_sorter;
  typedef typename CGAL::Edge_sorter<Nef_polyhedron, std::greater<Point_3>,
                                     std::greater<FT>, Container> Edge_sorter2;

  //  CGAL_NEF_SETDTHREAD(47*229*233*239);
  //  CGAL_NEF_SETDTHREAD(223);

  External_structure_builder esb;
  SFace_separator sf_sep;
  N.delegate(sf_sep,false, false);

  Reflex_edge_searcher res(Sphere_point(1,0,0));
  N.delegate(res,false,false);
  
  std::cerr << "number of reflex sedges"
	    << std::distance(res.positive_redges_begin(), 
			     res.positive_redges_end())
	    << ","
	    << std::distance(res.negative_redges_begin(), 
			     res.negative_redges_end())
	    << std::endl;

    Edge_sorter es(res.get_negative_redges());
    N.delegate(es);

    //    CGAL_NEF_SETDTHREAD(229*239*331);

    Reflex_sedge_iterator rei;
    for(rei=res.negative_redges_begin(); rei!=res.negative_redges_end(); ++rei) {
      Halfedge_handle e = (*rei);

      /*
      std::cerr << "handle negative reflex edge "
		<< ": " << e->source()->point() 
		<< "->" << e->twin()->source()->point() << std::endl;
      */

      Single_wall W(e,Vector_3(-1,0,0));
      if(!W.need_to_create_wall()) continue;
      
      Reflex_vertex_searcher rvs(Sphere_point(1,0,0));
      //      if((rvs.is_reflex_vertex(e->source())&2) == 2) {
      if(rvs.need_to_shoot(e, true)) {
	Ray_hit2 rh2a(Vector_3(-1,0,0), e->source());
	N.delegate(rh2a);
      }
      //      if((rvs.is_reflex_vertex(e->twin()->source())&2) == 2) {
      if(rvs.need_to_shoot(e->twin(), true)) {
	Ray_hit2 rh2a(Vector_3(-1,0,0), e->twin()->source());
	N.delegate(rh2a);
      }  
    }

    int i=0;
    for(rei=res.negative_redges_begin(); rei!=res.negative_redges_end(); ++rei) {
      Halfedge_handle e = (*rei);
      if((++i%100)==0)
	std::cerr << "handle negative reflex edge " << i << std::endl;
	  //		  << ": " << e->source()->point() 
	  //		  << "->" << e->twin()->source()->point() << std::endl;
      if(e->point().hx() > 0)
	e = e->twin();
      Single_wall W(e,Vector_3(-1,0,0));
      if(!W.need_to_create_wall()) continue;    
      N.delegate(W);
    }
    //    CGAL_NEF_SETDTHREAD(229*233*239);
    N.delegate(esb);
    N.delegate(res, false, false);

    std::cerr << "number of reflex sedges" 
	      << std::distance(res.positive_redges_begin(), 
			       res.positive_redges_end())
	      << ","
	      << std::distance(res.negative_redges_begin(), 
			       res.negative_redges_end())
	      << std::endl;

    CGAL_assertion(N.is_valid(0,0));

    Reflex_edge_searcher& res2 = res;

    Edge_sorter2 es2(res2.get_positive_redges());
    N.delegate(es2);

    //    Reflex_sedge_iterator rei;
    for(rei=res2.positive_redges_begin(); rei!=res2.positive_redges_end(); ++rei) {
      Halfedge_handle e = (*rei);

      CGAL_assertion(e->source()->point() >
		     e->twin()->source()->point());
      Single_wall W(e,Vector_3(1,0,0));
      if(!W.need_to_create_wall()) continue;
      
      Reflex_vertex_searcher rvs(Sphere_point(1,0,0));
      //      if((rvs.is_reflex_vertex(e->source())&1) == 1) {
      if(rvs.need_to_shoot(e, false)) {
	//	  std::cerr << "shoot source" << std::endl;
	Ray_hit2 rh2a(Vector_3(1,0,0), e->source());
	N.delegate(rh2a);
      }
      //      if((rvs.is_reflex_vertex(e->twin()->source())&1) == 1) {
      if(rvs.need_to_shoot(e->twin(), false)) {
	//	  std::cerr << "shoot target" << std::endl;
	Ray_hit2 rh2a(Vector_3(1,0,0), e->twin()->source());
	N.delegate(rh2a);
      }
    }
    
    i=0;
    for(rei=res2.positive_redges_begin(); rei!=res2.positive_redges_end(); ++rei) {
      Halfedge_handle e = (*rei);
      if((++i%100)==0)
	std::cerr << "handle positive reflex edge " << i << std::endl;
	  //		  << ": " << e->source()->point() 
	  //		  << "->" << e->twin()->source()->point() << std::endl;
      Single_wall W(e,Vector_3(1,0,0));
      if(!W.need_to_create_wall()) continue;
      N.delegate(W);
    }
    
    N.delegate(esb);    
    CGAL_assertion(N.is_valid(0,0));
    //    N.delegate(res2, false, false);

    std::cerr << "number of reflex sedges" 
	      << std::distance(res2.positive_redges_begin(), 
			       res2.positive_redges_end())
	      << ","
	      << std::distance(res2.negative_redges_begin(), 
			       res2.negative_redges_end())
	      << std::endl;
   
  typename Nef_polyhedron::Halfedge_const_iterator ei;
  CGAL_forall_halfedges(ei, N) {
    if(ei->out_sedge() == typename Nef_polyhedron::SHalfedge_const_handle())
      std::cerr << "isolated edge " << ei->source()->point() 
		<< "->" << ei->twin()->source()->point() << std::endl;
  }

  //  CGAL_NEF_SETDTHREAD(43*227*229*233);

#ifndef CGAL_MINKOWSKI3_OLD_DECOMPOSITION
  YVertical_wall_builder2 Y;
  N.delegate(Y,false,false);
#else
  YVertical_wall_builder Y;
  N.delegate(Y,false,false);

  std::cerr << "number of yreflex sedges "
	    << std::distance(Y.pos_begin(), 
			     Y.pos_end())
	    << ","
	    << std::distance(Y.neg_begin(), 
			     Y.neg_end())
	    << std::endl;
  

  typename YVertical_wall_builder::Vertical_redge_iterator vri=Y.pos_begin();
  for(; vri != Y.pos_end(); ++vri) {
    Single_wall2 W((*vri),Sphere_point(0, 1,0));
    N.delegate(W);
  }

  for(vri=Y.neg_begin(); vri != Y.neg_end(); ++vri) {
    Single_wall2 W((*vri),Sphere_point(0,-1,0));
    N.delegate(W);
  }
#endif

  N.delegate(esb);

  CGAL_assertion_code(typename Nef_polyhedron::SHalfedge_const_iterator cse);
  CGAL_assertion_code(CGAL_forall_shalfedges(cse, N)
    if(cse->incident_sface()->mark())
      CGAL_assertion(!CGAL::is_reflex_sedge_in_any_direction<Nef_polyhedron>(cse)));
}
