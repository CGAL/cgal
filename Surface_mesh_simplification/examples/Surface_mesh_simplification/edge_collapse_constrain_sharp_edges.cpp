#include <iostream>
#include <fstream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Constrained_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Unique_hash_map.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <CGAL/property_map.h>
#include <cmath>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_3 Point_3;
typedef CGAL::Polyhedron_3<Kernel> Surface_mesh;
typedef boost::graph_traits<Surface_mesh const>::edge_descriptor edge_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification ;


typedef Surface_mesh::Facet_iterator Facet_iterator;
typedef Surface_mesh::Halfedge_handle Halfedge_handle;
typedef Surface_mesh::Halfedge_iterator Halfedge_iterator;


//
// BGL property map which indicates whether an edge is marked as non-removable
//
struct Constrained_edge_map : public boost::put_get_helper<bool,Constrained_edge_map>
{
  typedef boost::readable_property_map_tag      category;
  typedef bool                                  value_type;
  typedef bool                                  reference;
  typedef edge_descriptor                       key_type;

  Constrained_edge_map(const CGAL::Unique_hash_map<key_type,bool>& aConstraints)
    : mConstraints(aConstraints) {}

  reference operator[](key_type const& e) const { return  is_constrained(e); }

  bool is_constrained( key_type const& e ) const {
    return mConstraints.is_defined(e) ? mConstraints[e] : false ; }

private:
  const CGAL::Unique_hash_map<key_type,bool>& mConstraints;
};

int main( int argc, char** argv )
{
  CGAL::Unique_hash_map<edge_descriptor,bool> constraint_hmap(false);

  Surface_mesh surface_mesh;

  if (argc!=2){
    std::cerr<< "Usage: " << argv[0] << " input.off\n";
    return 1;
  }

  std::ifstream is(argv[1]);
  if(!is){
    std::cerr<< "Filename provided is invalid\n";
    return 1;
  }

  is >> surface_mesh  ;

  Constrained_edge_map constraints_map(constraint_hmap);
  SMS::Constrained_placement<SMS::Midpoint_placement<Surface_mesh>,
                             Constrained_edge_map > placement(constraints_map);

  // map used to check that constrained_edges and the points of its vertices
  // are preserved at the end of the simplification
  // Warning: the computation of the diedral angle is only an approximation and can
  //          be far from the real value and could influence the detection of sharp
  //          edges after the simplification
  std::map<Surface_mesh::Halfedge_handle,std::pair<Point_3, Point_3> >constrained_edges;
  std::size_t nb_sharp_edges=0;

  // detect sharp edges
  std::ofstream cst_output("constrained_edges.cgal");
  for(Surface_mesh::Edge_iterator eb = surface_mesh.edges_begin(), ee = surface_mesh.edges_end() ; eb != ee ; ++eb )
  {
    if ( eb->is_border_edge() ){
      ++nb_sharp_edges;
      constraint_hmap[eb]=true;
      constraint_hmap[eb->opposite()]=true;
      constrained_edges[eb]=std::make_pair( eb->opposite()->vertex()->point(),
                                            eb->vertex()->point() );
    }
    else{
      double angle = CGAL::Mesh_3::dihedral_angle(
        eb->opposite()->vertex()->point(),
        eb->vertex()->point(),
        eb->next()->vertex()->point(),
        eb->opposite()->next()->vertex()->point() );
      if ( std::abs(angle)<100 ){
        ++nb_sharp_edges;
        constraint_hmap[eb]=true;
        constraint_hmap[eb->opposite()]=true;
        constrained_edges[eb]=std::make_pair( eb->opposite()->vertex()->point(),
                                              eb->vertex()->point() );
        cst_output << "2 " << eb->opposite()->vertex()->point() << " "
                           << " "  << eb->vertex()->point() << "\n";
      }
    }
  }
  cst_output.close();

  // Contract the surface mesh as much as possible
  SMS::Count_stop_predicate<Surface_mesh> stop(0);

  int r
  = SMS::edge_collapse(surface_mesh
                       ,stop
                       ,CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index,surface_mesh))
                       .edge_index_map  (boost::get(CGAL::edge_external_index  ,surface_mesh))
                       .edge_is_constrained_map(constraints_map)
                       .get_placement(placement)
   );

  std::cout << "\nFinished...\n" << r << " edges removed.\n"
  << (surface_mesh.size_of_halfedges()/2) << " final edges.\n" ;
  std::ofstream os(argc > 2 ? argv[2] : "out.off") ; os << surface_mesh ;

  std::cout  << "Checking sharped edges were preserved...\n";
  // check sharp edges were preserved
  for(Surface_mesh::Edge_iterator eb = surface_mesh.edges_begin(), ee = surface_mesh.edges_end() ; eb != ee ; ++eb )
  {
    if ( eb->is_border_edge() ){
      --nb_sharp_edges;
      assert(
        constrained_edges[eb]==std::make_pair( eb->opposite()->vertex()->point(),
                                               eb->vertex()->point() ) );
    }
    else{
      double angle = CGAL::Mesh_3::dihedral_angle(
        eb->opposite()->vertex()->point(),
        eb->vertex()->point(),
        eb->next()->vertex()->point(),
        eb->opposite()->next()->vertex()->point() );
      if ( std::abs(angle)<100 ){
        ++nb_sharp_edges;
      assert(
        constrained_edges[eb]==std::make_pair( eb->opposite()->vertex()->point(),
                                               eb->vertex()->point() ) );
      }
    }
  }
  std::cout  << "OK\n";

  std::cout << "Check that no removable edge has been forgotten..." << std::endl;
  r = SMS::edge_collapse(surface_mesh
                         ,stop
                         ,CGAL::vertex_index_map(boost::get(CGAL::vertex_external_index, surface_mesh))
                         .edge_index_map  (boost::get(CGAL::edge_external_index, surface_mesh))
                         .edge_is_constrained_map(constraints_map)
                         .get_placement(placement)
   );

  assert(r==0);

  if ( r==0 )
    std::cout  << "OK\n";
  else{
    std::cout  << "ERROR! " << r << " edges removed!\n";
    return 1;
  }

  return 0;
}
