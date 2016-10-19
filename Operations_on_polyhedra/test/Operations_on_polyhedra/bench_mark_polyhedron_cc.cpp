#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Timer.h>
#include <CGAL/property_map.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

#include <CGAL/intersection_of_Polyhedra_3.h>
#include <CGAL/intersection_of_Polyhedra_3_refinement_visitor.h>
#include <CGAL/Polygon_mesh_processing/connected_components.h>

template <class Refs>
struct Halfedge_with_mark : public CGAL::HalfedgeDS_halfedge_base<Refs> {
  Halfedge_with_mark()
    : CGAL::HalfedgeDS_halfedge_base<Refs>(), m(false)
  {}

  bool m;   // A boundary marker for faces with different status
  void set_mark()
  {
    m = true;
  }
};

// An items type using my halfedge.
struct Items_with_mark_on_hedge : public CGAL::Polyhedron_items_3
{
  template <class Refs, class Traits>
  struct Halfedge_wrapper {
      typedef Halfedge_with_mark<Refs> Halfedge;
  };
};

template<class Polyhedron>
struct Edge_mark_property_map{
  typedef bool value_type;
  typedef value_type reference;
  typedef std::pair<typename Polyhedron::Halfedge_handle,Polyhedron*> key_type;
  typedef boost::read_write_property_map_tag category;

  friend reference get(Edge_mark_property_map,const key_type& key)
  {return key.first->m;}
  friend void put(Edge_mark_property_map,key_type key,value_type v)
  {key.first->m=v;}
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef CGAL::Polyhedron_3<Kernel,Items_with_mark_on_hedge> Polyhedron;

typedef CGAL::Node_visitor_refine_polyhedra<Polyhedron,CGAL::Default,
                                            Kernel,
                                            Edge_mark_property_map<Polyhedron>
  > Split_visitor;
typedef std::vector<Kernel::Point_3> Polyline_3;

struct Is_not_marked{
  bool operator()(Polyhedron::Halfedge_const_handle h) const{
    return !h->m;
  }
};

struct My_facet_marker{
  std::size_t pi;
  My_facet_marker():pi((std::numeric_limits<std::size_t>::max)()){}
  std::map<Polyhedron::Facet_const_handle, std::size_t> patch_ids;
  template <class Facet_handle_iterator>
  void mark(Facet_handle_iterator begin, Facet_handle_iterator end)
  {
    for (Facet_handle_iterator fit=begin; fit!=end; ++fit)
      patch_ids[*fit]=pi;
  }
  void start_new_connected_component(){++pi;}
};


std::size_t nb_cc_v1(Polyhedron& P, Is_not_marked adjacent, My_facet_marker& fm)
{
  CGAL::internal::corefinement::mark_connected_components(P, adjacent, fm);
  return fm.pi+1;
}


template <class Facet_id_pmap>
std::size_t nb_cc_v2(const Polyhedron& P, Is_not_marked adjacent,
                     Facet_id_pmap facet_ids)
{
  std::vector<std::size_t> patch_ids;
  std::vector<std::size_t> patch_sizes;
  CGAL::internal::corefinement::init_facet_indices(P, facet_ids);
  return CGAL::internal::corefinement::mark_connected_components_v2
    (P, adjacent, facet_ids, patch_ids, patch_sizes);
}

int main(int argc,char** argv) {

    if (argc!=3){
      std::cerr << "Usage "<< argv[0] << " file1.off file2.off\n";
      return 1;
    }

    Polyhedron P, Q;
    std::ifstream file(argv[1]);
    file >> P;
    file.close();
    file.open(argv[2]);
    file >> Q;
    file.close();

    CGAL::set_pretty_mode(std::cerr);

    std::cout << P.size_of_facets() << " "
              << Q.size_of_facets() << std::endl;

    std::list<Polyline_3> polylines;

    Split_visitor visitor;
    CGAL::Intersection_of_Polyhedra_3<Polyhedron,Kernel,Split_visitor>
      polyline_intersections(visitor);
    std::cout << "Vertices before " <<  P.size_of_vertices()
              << " " << Q.size_of_vertices() << std::endl;
    polyline_intersections( P,Q,std::back_inserter(polylines));
    std::cout << "Vertices after " <<  P.size_of_vertices()
              << " " << Q.size_of_vertices() << std::endl;

    std::cout << polylines.size() << " polylines" << std::endl;

    CGAL_assertion(P.is_valid());
    CGAL_assertion(Q.is_valid());

    Is_not_marked criterium;

    //using union-find
    My_facet_marker P_dfm;
    CGAL::Timer time;
    time.start();
    std::size_t P_nb_patches_v1=nb_cc_v1(P, criterium, P_dfm);
    time.stop();
    std::cout << "P v1 " << time.time() << " " << P_nb_patches_v1 << std::endl;
    time.reset();
    My_facet_marker Q_dfm;
    time.start();
    std::size_t Q_nb_patches_v1=nb_cc_v1(Q, criterium, Q_dfm);
    time.stop();
    std::cout << "Q v1 " << time.time() << " " << Q_nb_patches_v1 << std::endl;
    time.reset();

    //using propagation
    typedef std::map<Polyhedron::Facet_const_handle,std::size_t> Patch_id_map;
    Patch_id_map P_patch_map;
    time.start();
    std::size_t P_nb_patches_v2 =
      nb_cc_v2( P,
                criterium,
                boost::associative_property_map<Patch_id_map>(P_patch_map) );
    time.stop();
    std::cout << "P v2 " << time.time() << " " << P_nb_patches_v2 << std::endl;
    time.reset();

    Patch_id_map Q_patch_map;
    time.start();
    std::size_t Q_nb_patches_v2 =
      nb_cc_v2(Q, criterium,
               boost::associative_property_map<Patch_id_map>(Q_patch_map));
    time.stop();
    std::cout << "Q v2 " << time.time() << " " << Q_nb_patches_v2 << std::endl;
    time.reset();

    return 0;
}
