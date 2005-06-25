// file: examples/Voronoi_diagram_adaptor_2/connected_components.C
#include <CGAL/basic.h>

// standard includes
#include <iostream>
#include <fstream>
#include <cassert>

// includes for the Apollonius graph and the Voronoi diagram adaptor
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_filtered_traits_2.h>
#include <CGAL/Voronoi_diagram_adaptor_2.h>
#include <CGAL/Apollonius_graph_Voronoi_traits_2.h>

// typedefs for defining the adaptor
typedef CGAL::Simple_cartesian<double>                K;
typedef CGAL::Apollonius_graph_filtered_traits_2<K>   Gt;
typedef CGAL::Apollonius_graph_2<Gt>                  AG;
typedef CGAL::Apollonius_graph_Voronoi_traits_2<AG>   VT;
typedef CGAL::Voronoi_diagram_adaptor_2<AG,VT>        VD;

// a generic functor that computes the connected components
template<class VDA_t>
class Connected_components
{
 private:
  typedef VDA_t                                            VDA;
  typedef typename VDA::Halfedge_iterator                  Halfedge_iterator;
  typedef typename VDA::Halfedge_handle                    Halfedge_handle;
  typedef typename VDA::Vertex_handle                      Vertex_handle;
  typedef typename VDA::Halfedge                           Halfedge;
  typedef typename VDA::Vertex                             Vertex;
  typedef typename VDA::Halfedge_around_vertex_circulator  HAVC;

 public:
  typedef VDA  Voronoi_diagram_adaptor_2;

  typedef typename Voronoi_diagram_adaptor_2::size_type  size_type;
  typedef size_type                                      result_type;
  typedef Voronoi_diagram_adaptor_2                      argument_type;

  struct Arity { enum { arity = 1 }; };

 private:
  struct Halfedge_less {
    bool operator()(const Halfedge& h1,	const Halfedge& h2) const {
      typename Voronoi_diagram_adaptor_2::Dual_edge e1 = h1.dual_edge();
      typename Voronoi_diagram_adaptor_2::Dual_edge e2 = h2.dual_edge();

      if ( e1.first != e2.first ) { return e1.first < e2.first; }
      return e1.second < e2.second;
    }
  };

  typedef std::map<Halfedge,bool,Halfedge_less>  Halfedge_map;

  void mark_edge(const Halfedge& e, Halfedge_map& e_map) const {
    e_map[e] = true;
    e_map[*e.opposite()] = true;
  }

  bool is_unmarked(const Halfedge& e, const Halfedge_map& e_map) const {
    return e_map.find(e) == e_map.end();
  }

  void dfs(const Voronoi_diagram_adaptor_2& vda, const Halfedge& e,
	      Halfedge_map& e_map) const
  {
    Halfedge e_opp = *e.opposite();
    mark_edge(e, e_map);

    HAVC ec, ec_start;
    if ( e.has_source() ) {
      ec = vda.incident_halfedges(e.source());
      ec_start = ec;
      do {
        if ( *ec != e && *ec != e_opp && is_unmarked(*ec, e_map) ) {
          dfs(vda, *ec, e_map);
        }
        ec++;
      } while (ec != ec_start);
    }

    if ( e.has_target() ) {
      ec = vda.incident_halfedges(e.target());
      ec_start = ec;
      do {
        if ( *ec != e && *ec != e_opp && is_unmarked(*ec, e_map) ) {
          dfs(vda, *ec, e_map);
        }
        ec++;
      } while (ec != ec_start);
    }
  }

 public:
  size_type operator()(const Voronoi_diagram_adaptor_2& vda) const
  {
    Halfedge_map e_map;

    size_type n_components = 0;
    for (Halfedge_iterator eit = vda.halfedges_begin();
         eit != vda.halfedges_end(); ++eit) {
      if ( is_unmarked(*eit, e_map) ) {
        n_components++;
        dfs(vda, *eit, e_map);
      }
    }
    return n_components;
  }
};


void eval_num_of_connected_components(char* fname) {
  std::ifstream ifs(fname);
  assert( ifs );

  std::cout << std::endl << "Input sites:" << std::endl;

  AG ag;
  AG::Site_2 s;
  while ( ifs >> s ) {
    std::cout << s << std::endl;
    ag.insert(s);
  }
  ifs.close();

  VD vd(ag);

  std::cout << std::endl << "# of connected components: "
	    << Connected_components<VD>()(vd) << std::endl << std::endl;

  // compute the number of bisectong segments, bisecting rays and
  // full bisectors in the Voronoi diagram
  VD::Edge_iterator eit;
  VD::size_type n_bisectors = 0, n_rays = 0, n_segments = 0;
  for (eit = vd.edges_begin(); eit != vd.edges_end(); ++eit) {
    if ( eit->curve().is_bisector() ) { n_bisectors++; }
    else if ( eit->curve().is_ray() ) { n_rays++; }
    else { n_segments++; }
  }

  std::cout << "# of bisecting segments: " << n_segments << std::endl;
  std::cout << "# of bisecting rays:     " << n_rays << std::endl;
  std::cout << "# of full bisectors:     " << n_bisectors << std::endl;
  std::cout << std::endl;
  std::cout << "===========================================" << std::endl;
}


int main()
{
  eval_num_of_connected_components("data/data1.ag.cin");
  eval_num_of_connected_components("data/data2.ag.cin");
  eval_num_of_connected_components("data/data3.ag.cin");
  eval_num_of_connected_components("data/degenerate.ag.cin");
  eval_num_of_connected_components("data/1D.ag.cin");

  return 0;
}
