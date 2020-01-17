//#define POLY 

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/boost/graph/helpers.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <CGAL/assertions.h>

#include <CGAL/boost/graph/Euler_operations.h>

#include <cassert>
#include <vector>
#include <set>
#include <fstream>


#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
namespace CGAL{
namespace Polygon_mesh_processing{

/*!
\brief triangulate_hole_with_cdt_2

\tparam PolygonMesh a model of `MutableFaceGraph`
\tparam NamedParameters a sequence of \ref pmp_namedparameters "Named Parameters"
\param pmesh polygon mesh which has the hole
\param border_halfedge a border halfedge incident to the hole
\param np optional sequence of \ref pmp_namedparameters "Named Parameters" among the ones listed below

\cgalNamedParamsBegin
   \cgalParamBegin{vertex_point_map} the property map with the points associated to the vertices of `pmesh`.
       If this parameter is omitted, an internal property map for
       `CGAL::vertex_point_t` should be available in `PolygonMesh`
       \cgalParamEnd
   \cgalParamBegin{geom_traits} a geometric traits class instance \cgalParamEnd
\cgalNamedParamsEnd

* \return `true` if the hole has been triangulated, `false` otherwise.
*/
template<typename PolygonMesh,
         typename NamedParameters>
bool
triangulate_hole_with_cdt_2(PolygonMesh& pmesh,
                            typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
                            const NamedParameters& np)
{
  typedef Halfedge_around_face_circulator<PolygonMesh>   Hedge_around_face_circulator;
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;
  typedef typename GetVertexPointMap<PolygonMesh,
      NamedParameters>::const_type Vpm;
  typedef typename Kernel_traits<
      typename boost::property_traits<Vpm>::value_type >::Kernel Kernel;
  typedef typename Kernel::Point_3 Point_3;
  typedef typename Kernel::Point_2 Point_2;
  typedef typename boost::graph_traits<PolygonMesh>::face_descriptor face_descriptor;
  Vpm vpm = parameters::choose_parameter(parameters::get_parameter(np, internal_np::vertex_point),
                                         get_const_property_map(boost::vertex_point, pmesh));
  
  std::vector<Point_3> pts;
  
  Hedge_around_face_circulator circ(border_halfedge,pmesh), done(circ);
  face_descriptor new_face = add_face(pmesh);
  set_halfedge(new_face, border_halfedge, pmesh);
  
  do{
    pts.push_back(get(vpm, target(*circ, pmesh)));
    set_face(*circ, new_face, pmesh);
  } while (++circ != done);
  
  typename Kernel::Plane_3 plane;
  linear_least_squares_fitting_3(pts.begin(),pts.end(),plane,CGAL::Dimension_tag<0>());
  
  //Project on plane
  std::vector<Point_2> polyline_2d;
  polyline_2d.reserve(pts.size());
  
  for(const auto& p : pts)
  {
    polyline_2d.push_back(plane.to_2d(p));
  }
  if(!CGAL::is_simple_2(polyline_2d.begin(), polyline_2d.end(), Kernel()))
    return false;
  
  return triangulate_face(new_face, pmesh, np);
}

}}
typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic;


template <class Polyhedron>
void read_poly(const char* file_name, Polyhedron& poly) {
  poly.clear();
  
  std::ifstream input(file_name);
  if ( !input || !(input >> poly)  || (num_vertices(poly) == 0)){
    std::cerr << "  Error: can not read file." << std::endl;
    assert(false);
  }
}

template <class Polyhedron, class Halfedge_handle>
void detect_borders(Polyhedron& poly, std::vector<Halfedge_handle>& border_reps)
{
  typedef CGAL::Halfedge_around_face_circulator<Polyhedron> Halfedge_around_facet_circulator;
  border_reps.clear();
  std::set<Halfedge_handle> border_map;
  for(Halfedge_handle h :  halfedges(poly)){
    if(face(h,poly)== boost::graph_traits<Polyhedron>::null_face() && border_map.find(h) == border_map.end()){
      border_reps.push_back(h);
      Halfedge_around_facet_circulator hf_around_facet(h,poly), done(hf_around_facet);
      do {
        bool insertion_ok = border_map.insert(*hf_around_facet).second;
        assert(insertion_ok);
      } while(++hf_around_facet != done);
    }
  }
}

template <class Polyhedron, class Halfedge_handle>
void read_poly_with_borders(const char* file_name, Polyhedron& poly, std::vector<Halfedge_handle>& border_reps) 
{
  read_poly(file_name, poly);
  detect_borders(poly, border_reps);
}

/******************************************************************/
template <class Polyhedron>
void test_triangulate_hole_with_cdt_2(const char* file_name) {
  typedef typename boost::graph_traits<Polyhedron>::face_descriptor Facet_handle;
  typedef typename boost::graph_traits<Polyhedron>::halfedge_descriptor Halfedge_handle;
  
  std::cout << "test_triangulate_hole:" << std::endl;
  std::cout << "  File: "<< file_name  << std::endl;
  Polyhedron poly;
  std::vector<Halfedge_handle> border_reps;
  read_poly_with_borders(file_name, poly, border_reps);
  
  for(typename std::vector<Halfedge_handle>::iterator it = border_reps.begin(); it != border_reps.end(); ++it) {
    std::vector<Facet_handle> patch;  
    if(!CGAL::Polygon_mesh_processing::triangulate_hole_with_cdt_2(poly, *it, CGAL::parameters::output_iterator(std::back_inserter(patch))))
    {
      std::cerr << " Error: border is not simple." << std::endl;
      assert(false);
    }
    if(patch.empty())
    {
      std::cerr << " Error: No patch generated." << std::endl;
      assert(false);
    }
  }
  
  if(!poly.is_valid() || ! is_closed(poly)) {
    std::cerr << "  Error: patched polyhedron is not valid or closed." << std::endl;
    assert(false);
  }
  
  std::ofstream out("out.off");
  out << poly;
  out.close();         
  std::cout << "  Done!" << std::endl;
}

template <class Kernel>
void test_hole_filling() {
  typedef CGAL::Surface_mesh<typename Kernel::Point_3> Polyhedron;
  
  std::vector<std::string> input_files;
  input_files.push_back("data/w.off");
  //for(std::vector<std::string>::iterator it = input_files.begin(); it != input_files.end(); ++it) {
    test_triangulate_hole_with_cdt_2<Polyhedron>("/home/gimeno/Data/tmp/triforce.off");
          //it->c_str());
  //}
}

int main() 
{
  test_hole_filling<Epic>();
  std::cout << "All Done!" << std::endl;
}
