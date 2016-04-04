#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/HalfedgeDS_face_max_base_with_id.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Timer.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <CGAL/utility.h>

#include <CGAL/intersection_of_Polyhedra_3.h>
#include <CGAL/intersection_of_Polyhedra_3_refinement_visitor.h>
#include <CGAL/internal/corefinement/Polyhedra_output_builder.h>
#include <CGAL/iterator.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/selection.h>

template <class Refs>
struct Halfedge_with_mark : public CGAL::HalfedgeDS_halfedge_base<Refs> {
  Halfedge_with_mark()
    : CGAL::HalfedgeDS_halfedge_base<Refs>(), m(false)
  {}

  bool m;   // A boundary marker for faces with different status
  void set_mark(bool v)
  {
    m = v;
  }
};

template<class Polyhedron>
struct Edge_mark_property_map4coref{
  typedef bool value_type;
  typedef value_type reference;
  typedef std::pair<typename Polyhedron::Halfedge_handle,Polyhedron*> key_type;
  typedef boost::read_write_property_map_tag category;

  friend reference get(Edge_mark_property_map4coref,const key_type& key) {return key.first->m;}
  friend void put(Edge_mark_property_map4coref,key_type key,value_type v) {key.first->m=v;}
};

template<class Polyhedron>
struct Edge_mark_property_map{
  typedef bool value_type;
  typedef value_type reference;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor key_type;
  typedef boost::readable_property_map_tag category;

  friend reference get(Edge_mark_property_map,const key_type& key) {return key.halfedge()->m;}
};

// An items type using my halfedge.
struct Items_with_mark_on_hedge : public CGAL::Polyhedron_items_3
{
  template <class Refs, class Traits>
  struct Halfedge_wrapper {
      typedef Halfedge_with_mark<Refs> Halfedge;
  };

  template < class Refs, class Traits>
  struct Face_wrapper {
      typedef CGAL::HalfedgeDS_face_max_base_with_id< Refs, CGAL::Tag_false, std::size_t>  Face;
  };
  //~ template < class Refs, class Traits>
  //~ struct Vertex_wrapper {
      //~ typedef typename Traits::Point_3 Point;
      //~ typedef CGAL::HalfedgeDS_vertex_max_base_with_id< Refs, Point, std::size_t> Vertex;
  //~ };
};



typedef CGAL::Exact_predicates_inexact_constructions_kernel              Kernel;
typedef CGAL::Polyhedron_3<Kernel,Items_with_mark_on_hedge>          Polyhedron;
typedef Edge_mark_property_map4coref<Polyhedron> Edge_ppmap;
typedef std::map<Polyhedron::Facet_const_handle,std::size_t>       Facet_id_map;
typedef boost::associative_property_map<Facet_id_map>             Facet_id_pmap;
typedef CGAL::Corefinement
            ::Polyhedra_output_builder< Polyhedron,
                                        Facet_id_pmap,
                                        CGAL::Default,
                                        CGAL::Default,
                                        Edge_ppmap >             Output_builder;
typedef CGAL::Node_visitor_refine_polyhedra<Polyhedron,
                                            Output_builder,
                                            Kernel,
                                            Edge_ppmap>           Split_visitor;


typedef Polyhedron::Face_handle                                     Face_handle;

struct Is_facet_selected{
  std::set<Face_handle>& selection;

  Is_facet_selected(std::set<Face_handle>& sel)
    :selection(sel)
  {}

  friend bool get(Is_facet_selected is, Face_handle f){
    return is.selection.count(f);
  }

  friend void put(Is_facet_selected is, Face_handle f, bool b){
    if (b)
      is.selection.insert(f);
    else
      is.selection.erase(f);
  }
};

namespace PMP=CGAL::Polygon_mesh_processing;

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

  CGAL::Emptyset_iterator output_it;
  Facet_id_map P_facet_id_map, Q_facet_id_map;

  CGAL::cpp11::array<boost::optional<Polyhedron*>, 4 > desired_output;
  Polyhedron inter, union_;
  desired_output[Output_builder::P_UNION_Q]=boost::make_optional( &union_ );

  Output_builder output_builder(P, Q,
                                desired_output,
                                Facet_id_pmap(P_facet_id_map),
                                Facet_id_pmap(Q_facet_id_map) );
  Split_visitor visitor(output_builder);

  CGAL::Intersection_of_Polyhedra_3<Polyhedron,Kernel,Split_visitor> polyline_intersections(visitor);
  std::cout << "Vertices before " <<  P.size_of_vertices()
            << " " << Q.size_of_vertices() << std::endl;
  polyline_intersections(P, Q, output_it);
  std::cout << "Vertices after " <<  P.size_of_vertices()
            << " " << Q.size_of_vertices() << std::endl;

  if ( output_builder.union_valid() ) std::cout << "P-Q is valid\n";
  else{
    std::cout << "P-Q is invalid\n";
    return 1;
  }

  CGAL_assertion(union_.is_valid());

  std::ofstream output("P_union_Q.off");
  output << std::setprecision(17) << union_;
  output.close();

  std::size_t i=0;
  BOOST_FOREACH(Polyhedron::Face_handle f, faces(union_))
    f->id()=i++;

  std::set<Polyhedron::Face_handle> faces;
  BOOST_FOREACH(Polyhedron::Halfedge_handle h, halfedges(union_))
    if( !h->is_border() && h->m){
      // insert all faces incident to the target vertex
      BOOST_FOREACH(Polyhedron::Halfedge_handle h2, halfedges_around_target(h,union_))
        if (!h2->is_border())
          faces.insert(h2->face());
    }

  std::cout << "faces.size() " << faces.size() << "\n";
  CGAL::expand_face_selection(faces, P, 2, Is_facet_selected(faces), CGAL::Emptyset_iterator());

  {
    std::ofstream output("faces.selection.txt");
    output << "\n";
    BOOST_FOREACH(Polyhedron::Face_handle fh, faces)
      output << fh->id() << " ";
    output << "\n";

    std::set<std::pair<Kernel::Point_3,Kernel::Point_3> > csts;
    BOOST_FOREACH(boost::graph_traits<Polyhedron>::edge_descriptor ed, edges(union_))
    {
      if (ed.halfedge()->m) csts.insert(CGAL::make_sorted_pair(ed.halfedge()->vertex()->point(), ed.halfedge()->opposite()->vertex()->point()));
    }

    std::ofstream tmp("/tmp/tmp.off");
    tmp << std::setprecision(17);
    tmp << union_;
    tmp.close();
    std::ifstream in("/tmp/tmp.off");
    Polyhedron Pp;
    in >> Pp;

    int i=0;
    BOOST_FOREACH(boost::graph_traits<Polyhedron>::edge_descriptor ed, edges(Pp))
    {
      if (csts.count(CGAL::make_sorted_pair(ed.halfedge()->vertex()->point(), ed.halfedge()->opposite()->vertex()->point())))
        output << i << " ";
      ++i;
    }
    output << "\n";
  }

  //~ CGAL::Polygon_mesh_processing::isotropic_remeshing(P, faces,0.1);
  CGAL::Polygon_mesh_processing::isotropic_remeshing(
    faces,0.3, union_,
    PMP::parameters::edge_is_constrained_map(Edge_mark_property_map<Polyhedron>()).number_of_iterations(15).smooth_along_features(true)
  );

  output.open("P_union_Q_remeshed.off");
  output << std::setprecision(17) << union_;


  return 0;
}
