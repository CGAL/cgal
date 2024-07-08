#pragma once

#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Side_of_triangle_mesh.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>
#include <CGAL/Aff_transformation_3.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Bbox_3.h>
#include <functional>
#include "utils.h"

#include "../query_replace/cmap_query_replace.h"
#include "../query_replace/init_to_preserve_for_query_replace.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT       FT;
typedef Kernel::Point_3  Point;
typedef Kernel::Vector_3 Vector;

typedef typename Kernel::Triangle_3 Triangle;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> Tree;

class DartInfo
{
public:
  struct FaceAttrValue { bool explored = false; };
  struct VolumeAttrValue { char template_id = 0; char iteration = 0; };
  template < class Refs >
  struct Dart_wrapper
  {
    typedef CGAL::Cell_attribute<Refs, FaceAttrValue> FaceAttr;
    typedef CGAL::Cell_attribute<Refs, VolumeAttrValue> VolumeAttr;
    typedef std::tuple<CGAL::Cell_attribute_with_point<Refs>, void, FaceAttr, VolumeAttr> Attributes;
  };
};
typedef CGAL::Linear_cell_complex_traits<3,Kernel> Mytraits;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3, Mytraits, DartInfo> LCC;
typedef typename LCC::Dart_handle Dart_handle;
typedef typename LCC::Vertex_attribute_handle Vertex_handle;
typedef typename LCC::size_type   size_type;


const size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

// to be removed
size_type debug_node_mark;
size_type debug_edge_mark;

namespace CGAL::HexRefinement::TwoRefinement {
  enum Plane { XY, YZ, ZX, YX = XY, ZY = YZ, XZ = ZX };
  struct Grid { 
    Point center;
    double size;
    int nb_subdiv_per_dim;
  };

  struct HexMeshingData { 
    LCC lcc;
    size_type identified_mark, template_mark, corner_mark;
    Pattern_substituer<LCC> regular_templates, partial_templates;
    Grid grid;

    std::array<std::vector<Dart_handle>, 3> first_face_of_planes;
  };

  struct PlaneData {
    Plane plane;

    std::vector<Dart_handle> marked_nodes;
    std::vector<Dart_handle> faces_of_planes;
    std::vector<Dart_handle> volumes_to_refine;
    std::vector<Dart_handle> partial_templates_to_refine;
  };

  template <unsigned int i>
  auto get_or_create_attr(LCC& lcc, Dart_handle dart){
    auto attr = lcc.attribute<i>(dart);
    
    if (attr == nullptr){
      attr = lcc.create_attribute<i>();
      lcc.set_attribute<i>(dart, attr);
    }
    return attr;
  }

  Dart_handle find_3_template_origin(LCC& lcc, Dart_handle marked_face, size_type corner_mark) {

    Dart_handle dart = marked_face;

    // Get the origin dart : Find the two unmarked node on the face 
    // since the 3 template is created by adjacent two 2-templates
    bool found = false;
    for (int i = 0; i < 6; i++)
    {
      Dart_handle other_d_nonmark = lcc.beta(dart, 1, 1);
      if (!lcc.is_marked(dart, corner_mark) && !lcc.is_marked(other_d_nonmark, corner_mark))
      {
        found = true;
        break;
      }

      dart = lcc.beta(dart, 1);
    }

    assert(found);

    return dart;
  }

  template <unsigned int i>
  void mark_all_0_cells(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<i, 0>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++)
      if (!lcc.is_marked(dit, mark)) lcc.mark_cell<0>(dit, mark);
  }

  void clean_up_3_template(HexMeshingData &hdata, const Dart_handle &origin_dart, const Dart_handle &upper_edge, const Dart_handle lower_edge, Dart_handle &face1, Dart_handle &face2)
  {
    LCC& lcc = hdata.lcc;

    Dart_handle lower_mid_1 = lcc.insert_barycenter_in_cell<1>(lower_edge);
    Dart_handle lower_mid_2 = lcc.beta(lower_mid_1, 2, 1);

    // Contract the lower edge asweel, this leads to two faces being partially contracted.
    lcc.contract_cell<1>(lower_mid_1);
    lcc.contract_cell<1>(lower_mid_2);

    assert(!lcc.is_marked(lower_mid_2, hdata.corner_mark));

    // Contract the two remaining 2-darts faces
    Dart_handle face_to_remove_1 = lcc.beta(upper_edge, 2, 1);
    Dart_handle face_to_remove_2 = lcc.beta(upper_edge, 2, 1, 1, 2, 1, 2);

    assert(lcc.darts_of_orbit<1>(face_to_remove_1).size() == 2);
    assert(lcc.darts_of_orbit<1>(face_to_remove_2).size() == 2);

    // Contract the two remaining 2-darts faces
    lcc.contract_cell<2>(face_to_remove_1);
    lcc.contract_cell<2>(face_to_remove_2);

    // Remove the created volume from the partial query replace and sew the two neighboring volumes
    lcc.remove_cell<3>(origin_dart);

    if (face1 != lcc.null_dart_descriptor && face2 != lcc.null_dart_descriptor)
      lcc.sew<3>(face1, face2);
  }

  void refine_3_template(HexMeshingData &hdata, Dart_handle marked_face)
  {
    // TODO Might be written better 
    LCC& lcc = hdata.lcc;

    Dart_handle origin_dart = find_3_template_origin(lcc, marked_face, hdata.corner_mark);  
    Dart_handle vol2_origin_dart = lcc.beta(origin_dart, 3);

    Dart_handle upper_d1 = origin_dart;
    Dart_handle upper_d2 = lcc.beta(origin_dart, 1, 1);

    assert(!lcc.is_marked(upper_d2, hdata.corner_mark));

    Dart_handle upper_edge = lcc.insert_cell_1_in_cell_2(upper_d1, upper_d2);
    Dart_handle vol2_upper_edge = lcc.beta(upper_edge, 3); 

    // Query replace with the partial 3-template, making it into two volumes
    size_type p = hdata.partial_templates.query_replace_one_volume(lcc, origin_dart, hdata.template_mark);
    assert(p == 0);

    // Also replace the other connected volume that is 3 template
    p = hdata.partial_templates.query_replace_one_volume(lcc, lcc.beta(origin_dart, 3), hdata.template_mark);
    assert(p == 0);

    // Assert the origin dart have different directions on the adjacent volume
    // So assert that their opposite are equal
    assert(lcc.beta(origin_dart, 1) == lcc.beta(vol2_origin_dart, 0, 3));

    // Face of the two neighboring volumes to the created volume
    Dart_handle face1 = lcc.beta(origin_dart, 2, 3);
    Dart_handle face2 = lcc.beta(origin_dart, 1, 2, 3);
    Dart_handle lower_edge = lcc.beta(upper_edge, 2, 1, 1);

    Dart_handle vol2_face1 = lcc.beta(vol2_origin_dart, 2, 3);
    Dart_handle vol2_face2 = lcc.beta(vol2_origin_dart, 0, 2, 3); // 0 because opposite directions
    Dart_handle vol2_lower_edge = lcc.beta(vol2_upper_edge, 2, 1, 1);

    // Contract upper and lower edge into its barycenter 

    Dart_handle upper_mid_1 = lcc.insert_barycenter_in_cell<1>(upper_edge);
    Dart_handle upper_mid_2 = lcc.beta(upper_mid_1, 2, 1);
    
    assert(!lcc.is_marked(upper_mid_2, hdata.corner_mark));

    // contract the shared edge between the two volumes
    lcc.contract_cell<1>(upper_mid_1);
    lcc.contract_cell<1>(upper_mid_2);

    clean_up_3_template(hdata, origin_dart, upper_edge, lower_edge, face1, face2);
    clean_up_3_template(hdata, vol2_origin_dart, vol2_upper_edge, vol2_lower_edge, vol2_face1, vol2_face2);
  }

  void refine_marked_hexes(HexMeshingData& hdata, PlaneData& pdata)
  {
    LCC& lcc = hdata.lcc;
    int nbsub = 0; 
    int nb_3_tp = 0;

    for (auto& dart : pdata.faces_of_planes)
    {
      if (hdata.regular_templates.query_replace_one_face(lcc, dart, hdata.template_mark) != SIZE_T_MAX) nbsub++;
    }

    // Cannot easily assert if all faces has been correctly treated, because some faces don't have attr
    // and we don't refine 3/4 template faces.
    
    std::cout << nbsub << " face substitution was made" << std::endl;

    nbsub = 0;
    for (auto& dart : pdata.volumes_to_refine)
    {
      if (hdata.regular_templates.query_replace_one_volume(lcc, dart, hdata.template_mark) != SIZE_T_MAX) nbsub++;
    }

    // assert(nbsub == pdata.volumes_to_refine.size());

    // // Refine remaining 3 patterns
    for (auto marked_face : pdata.partial_templates_to_refine){
      refine_3_template(hdata, marked_face);
      nbsub++;
    }

    std::cout << nbsub << " volumic substitution was made" << std::endl;
  }


  Dart_handle __next_even_plane(LCC& lcc, Dart_handle start_plane, Plane plane){
    int f = plane == Plane::ZY ? 1 : 0;
    return lcc.beta(start_plane, f, 2, 3, 2, f, f, 2, 3, 2, f);
  }

  Dart_handle __first_vertex_of_even_plane(LCC& lcc, Plane plane) {
    switch (plane){
      case YZ: return lcc.beta(lcc.first_dart(), 0, 2, 1, 2, 3, 2, 1); 
      case XY: return lcc.beta(lcc.first_dart(), 2);
      case ZX: return lcc.beta(lcc.first_dart(), 0, 2, 0);
    }

    CGAL_assertion_msg(false, "Unexpected value for plane");
    return lcc.null_dart_descriptor;
  }

  /**
   * Mark 0-cells
   * Gather 2-cells adjacent to marked 0-cells
   * Gather 3-cells adjacents to marked 0-cells
   */
  void explore_face_of_plane(HexMeshingData& hdata, PlaneData& pdata, std::queue<Dart_handle>& queue, Dart_handle face, size_type explored_mark) {
    LCC& lcc = hdata.lcc;

    auto edges = lcc.darts_of_cell<2,1>(face);

    auto& face_attr = get_or_create_attr<2>(lcc, face)->info();

    if (face_attr.explored) return;

    // The 3-attr might not exist if the cell is not identified.
    int i = 0;
    
    // Add neighboring faces
    for (auto dit = edges.begin(), dend = edges.end(); dit != dend; dit++){
      bool explored = lcc.is_marked(dit, explored_mark);
      bool identified = lcc.is_marked(dit, hdata.identified_mark);
      bool has_opposite_face = !lcc.is_free<3>(lcc.beta(dit, 2));

      if (!explored){
        lcc.mark_cell<0>(dit, explored_mark);
      }

      if (!explored && identified){
        lcc.mark_cell<0>(dit, hdata.template_mark);
        pdata.marked_nodes.push_back(dit);
      }

      if (identified){
        // The 3attr is created only if one vertex is identiifed/marked
        get_or_create_attr<3>(lcc, face)->info().template_id++;
      }

      if (has_opposite_face) 
        queue.push(lcc.beta(dit, 2, 3, 2));

      ++i;
    }

    assert(i == 4);

    auto vol_handle = lcc.attribute<3>(face);
    face_attr.explored = true;
   
    if (vol_handle != nullptr) {
      assert(vol_handle->info().template_id > 0);
      pdata.faces_of_planes.push_back(face);
    }
  }

  void extract_darts_from_even_planes(HexMeshingData& hdata, PlaneData& pdata, Plane iterationPlane){
    LCC& lcc = hdata.lcc;
    pdata.plane = iterationPlane;

    size_type explored_mark = lcc.get_new_mark();
    size_type explored_faces = lcc.get_new_mark();
    
    for (auto start_plane : hdata.first_face_of_planes[iterationPlane]) {
      std::queue<Dart_handle> to_explore;
      to_explore.push(start_plane); // First face
      
      while (!to_explore.empty()) {
        Dart_handle front = to_explore.front();
        to_explore.pop();
        explore_face_of_plane(hdata, pdata, to_explore, front, explored_mark);
      }
    }
    
    // Could have also been done inside the loop
    std::vector<Dart_handle> additional_faces;

    for (auto face : pdata.faces_of_planes) {
      auto edges = lcc.darts_of_cell<2,0>(face);

      for (auto edge = edges.begin(); edge != edges.end(); edge++){
        // get incident faces to marked node to be refined
        // Incident faces normal to the plane 
        // We also need to prevent adding twice the faces by marking them
        
        if (!lcc.is_marked(edge, hdata.template_mark) && !lcc.is_marked(lcc.other_extremity(edge), hdata.template_mark))
          continue;

        auto top_face_1 = lcc.beta(edge, 2);
        auto top_face_2 = lcc.beta(edge, 3, 2);

        if (!lcc.is_whole_cell_marked<2>(top_face_1, explored_faces)){
          additional_faces.push_back(top_face_1);
          lcc.mark_cell<2>(top_face_1, explored_faces);
        }

        if (top_face_2 != lcc.null_dart_descriptor && !lcc.is_whole_cell_marked<2>(top_face_2, explored_faces)){
          additional_faces.push_back(top_face_2);
          lcc.mark_cell<2>(top_face_2, explored_faces);
        }
      }

      auto vol_handle = lcc.attribute<3>(face);

      // Also add the adjacent volumes if there is atleast one marked node
      // Because we are on odd/even layers, we can't accidently add twice a volume

      if (vol_handle == nullptr) continue; 

      auto& vol_attr = vol_handle->info();

      if (vol_attr.template_id != 3){
        pdata.volumes_to_refine.push_back(face);

        if (!lcc.is_free<3>(face)){
          Dart_handle other_face = lcc.beta(face, 3);
          auto& other_vol_attr = get_or_create_attr<3>(lcc, other_face)->info();
          other_vol_attr = vol_attr;
          pdata.volumes_to_refine.push_back(other_face);
        }
      }
      else {
        pdata.partial_templates_to_refine.push_back(face);
      }
    }

    pdata.faces_of_planes.insert(pdata.faces_of_planes.end(), additional_faces.begin(), additional_faces.end());
  }

  // Gather the faces first, if we don't we will have issues later manually 
  // Iterating to the next plane after a first refinement stage
  void extract_first_faces_of_planes(HexMeshingData& hdata){
    int nb_subdiv = hdata.grid.nb_subdiv_per_dim;

    for (int p = 0; p < 3; p++){
      Plane plane = (Plane)p;

      auto& start_of_planes = hdata.first_face_of_planes[p];
      
      Dart_handle start_plane = __first_vertex_of_even_plane(hdata.lcc, plane);
      
      for (int z = 0; z < nb_subdiv / 2; z++){
        start_of_planes.push_back(hdata.lcc.beta(start_plane, 0, 2));
        start_plane = __next_even_plane(hdata.lcc, start_plane, plane);
      }
    }
  }

  void mark_identified_cells_from_3_attrs(HexMeshingData& hdata) {
    LCC& lcc = hdata.lcc;

    auto identified = lcc.attributes<3>();

    for (auto it = identified.begin(), end = identified.end(); it != end; it++){
      mark_all_0_cells<3>(lcc, it->dart(), hdata.identified_mark);
    }
  }

  void propagation_stage(HexMeshingData& hdata, PlaneData& pdata){
    
  }

  void fix_adjacent_3_templates(HexMeshingData &hdata, Dart_handle &face)
  {
    LCC& lcc = hdata.lcc;

    auto &vol_attr = lcc.attribute<3>(face)->info();
    
    if (vol_attr.template_id != 3) return;

    auto edges = lcc.darts_of_cell<2, 1>(face);
    for (auto edge = edges.begin(); edge != edges.end(); edge++)
    {
      if (lcc.is_whole_cell_marked<1>(edge, hdata.template_mark))
        continue;

      // Transform nearby 3 templates into two 4 templates
      // if they both share the unmarked node

      auto other_vol_handle = lcc.attribute<3>(lcc.beta(edge, 2, 3, 2));

      if (other_vol_handle != nullptr && other_vol_handle->info().template_id != 3)
        continue;

      // Mark the unmarked node
      if (!lcc.is_marked(edge, hdata.template_mark))
        lcc.mark_cell<0>(edge, hdata.template_mark);

      auto other_ext = lcc.other_extremity(edge);
      if (!lcc.is_marked(other_ext, hdata.template_mark))
        lcc.mark_cell<0>(other_ext, hdata.template_mark);

      return;
    }
  }

  // If not all marks are chained consecutively, mark the whole face
  void fix_mark_connectivity(HexMeshingData &hdata, Dart_handle face){
    LCC& lcc = hdata.lcc;
    
    auto edges = lcc.darts_of_cell<2, 1>(face);
    bool connected = true;
    bool ended = false;
    bool started = false;

    for (auto edge = edges.begin(); edge != edges.end(); edge++) {
      if (ended && lcc.is_marked(edge, hdata.template_mark)) {
        connected = false;
        break;
      }

      if (started && !lcc.is_marked(edge, hdata.template_mark)) {
        ended = true;
      }

      if (!started && lcc.is_marked(edge, hdata.template_mark)){
        started = true;
      }
    }


    if (!connected) {

      for (auto edge = edges.begin(); edge != edges.end(); edge++) {
        if (lcc.is_marked(edge, hdata.template_mark)) continue;

        lcc.mark_cell<0>(edge, hdata.template_mark);

        // TODO
        // auto incident_faces = lcc.one_dart_per_incident_cell<0,2>(edge);
        // for (auto iface = incident_faces.begin(); iface != incident_faces.end(); iface++){
        //   auto &face_attr = lcc.attribute<2>(iface)->info();

        // }

      }
    }

  }

  void fix_impossible_cases(HexMeshingData &hdata, PlaneData &pdata)
  {
    for (auto face : pdata.partial_templates_to_refine){
      fix_adjacent_3_templates(hdata, face);
    }

    for (auto face : pdata.faces_of_planes){
      fix_mark_connectivity(hdata, face);
    }
  }



  Grid generate_grid(LCC& lcc, Point pos, double size, int nb_subdiv_per_dim) {
    double s = size / nb_subdiv_per_dim;

    // CGAL_precondition_msg((nb_subdiv_per_dim % 2) == 0, "nb_subdiv must be even");
    // or CGAL_precondition_msg((nb_subdiv_per_dim > 1) == 0, "grid should be atleast 2x2");

    int half_dim = nb_subdiv_per_dim / 2;

    for (int x = -half_dim; x < half_dim; x++) {
      for (int y = -half_dim; y < half_dim; y++) {
        for (int z = -half_dim; z < half_dim; z++) {

          double x1 = pos.x() + x * s, y1 = pos.y() + y * s, z1 = pos.z() + z * s;
          double x2 = pos.x() + (x+1)*s, y2 = pos.y() + (y+1)*s, z2 = pos.z() + (z+1)*s;

          lcc.make_hexahedron(Point(x1,y1,z1), Point(x2,y1,z1),
                              Point(x2,y2,z1), Point(x1,y2,z1),
                              Point(x1,y2,z2), Point(x1,y1,z2),
                              Point(x2,y1,z2), Point(x2,y2,z2));
        }
      }
    }

    lcc.sew3_same_facets();

    // Connected components must be equal to 1
    return Grid{.center = pos, .size = size, .nb_subdiv_per_dim = nb_subdiv_per_dim };
  }

  Grid generate_grid(LCC& lcc, Tree& aabb, int nb_subdiv_per_dim) {
    auto bbox = aabb.bbox(); 

    Point center = {bbox.xmin() + (bbox.x_span()/2), 
                    bbox.ymin() + (bbox.y_span()/2), 
                    bbox.zmin() + (bbox.z_span()/2)};

    double max_size = std::max(std::max(bbox.x_span(), bbox.y_span()), bbox.z_span()); 
    
    return generate_grid(lcc, center, max_size, nb_subdiv_per_dim);
  }

  void create_vertices_for_templates(HexMeshingData& hdata, PlaneData& pdata)
  {

    // 2 noeuds marqué l'un à coté de l'autre ne produit pas de sommet
    // 1 noeud marqué a coté d'un noeud non marqué produit un sommet

    // TODO A changer pour une itération sur les arretes ou pas selon comment l'algo va s'executer

    std::vector<Dart_handle> edges_to_subdivide;
    LCC& lcc = hdata.lcc;

    auto arrete_done = lcc.get_new_mark();

    for (auto dart : pdata.marked_nodes)
    {
      for (auto nit = lcc.one_dart_per_incident_cell<1, 0>(dart).begin(),
                nend = lcc.one_dart_per_incident_cell<1, 0>(dart).end();
          nit != nend;
          nit++)
      {
        if (lcc.is_marked(nit, arrete_done))
          continue;

        // If the node is next to an other marked node, we don't have to create vertices
        if (lcc.is_marked(lcc.beta<1>(nit), hdata.template_mark)){
          lcc.mark_cell<1>(nit, arrete_done);
          continue;
        }

        edges_to_subdivide.push_back(nit);
        lcc.mark_cell<1>(nit, arrete_done);
      }

      dart = lcc.beta<1>(dart);
    }

    for (Dart_handle dart : edges_to_subdivide)
    {
      lcc.insert_barycenter_in_cell<1>(dart);
    }

    lcc.free_mark(arrete_done);
  }

  void mark_1template_face(LCC& lcc, size_type mark){
    auto dart = lcc.one_dart_per_cell<3>().begin();
    lcc.mark_cell<0>(lcc.beta(dart, 0, 2, 0, 0), mark);
  }

  void mark_2template_face(LCC& lcc, size_type mark){
    auto dart = lcc.one_dart_per_cell<3>().begin();
    lcc.mark_cell<0>(dart, mark);
    lcc.mark_cell<0>(lcc.beta(dart, 1), mark);
  }

  void mark_1template_volume(LCC& lcc, size_type mark){
    auto first = lcc.first_dart();
    lcc.mark_cell<0>(lcc.beta(first, 1, 2, 3, 1, 2, 0), mark);
  }

  void mark_2template_volume(LCC& lcc, size_type mark){
    auto first = lcc.first_dart();
    
    auto node1 = lcc.beta(first, 0, 2, 3, 1, 1, 2, 0);
    auto node2 = lcc.beta(node1, 0);

    lcc.mark_cell<0>(node1, mark);
    lcc.mark_cell<0>(node2, mark);
  }

  void mark_3template_partial_volume(LCC& lcc, size_type mark){
    auto first = lcc.first_dart();

    auto m2 = lcc.get_new_mark();

    auto node1 = lcc.beta(first, 0, 0, 2, 1, 1, 2);
    auto node2 = lcc.beta(node1, 1);
    auto node3 = lcc.beta(node2, 1);

    lcc.mark_cell<0>(node1, mark);
    lcc.mark_cell<0>(node2, mark);
    lcc.mark_cell<0>(node3, mark);
  }

  void mark_4template_volume(LCC& lcc, size_type mark){
    auto first = lcc.first_dart();

    auto node1 = first;
    auto node2 = lcc.beta(node1, 0);
    auto node3 = lcc.beta(node2, 0);
    auto node4 = lcc.beta(node3, 0);

    lcc.mark_cell<0>(node1, mark);
    lcc.mark_cell<0>(node2, mark);
    lcc.mark_cell<0>(node3, mark);
    lcc.mark_cell<0>(node4, mark);
  }



  void load_patterns(Pattern_substituer<LCC> &regular_templates, Pattern_substituer<LCC>& partial_3_template) {
    // BUG Manque des opérateurs de déplacements pour Pattern, le reserve est un fix temporaire
    // Pour pouvoir charger les patterns correctement sans réallocation
    regular_templates.m_fpatterns.reserve(10);
    regular_templates.m_vpatterns.reserve(10);

    regular_templates.load_additional_fpattern(CGAL::data_file_path("hexmeshing/fpattern/pattern1-face.moka"), mark_1template_face);
    regular_templates.load_additional_fpattern(CGAL::data_file_path("hexmeshing/fpattern/pattern2-face.moka"), mark_2template_face);
    
    regular_templates.load_additional_vpattern(CGAL::data_file_path("hexmeshing/vpattern/complete/pattern1.moka"), mark_1template_volume);
    regular_templates.load_additional_vpattern(CGAL::data_file_path("hexmeshing/vpattern/complete/pattern2.moka"), mark_2template_volume);
    regular_templates.load_additional_vpattern(CGAL::data_file_path("hexmeshing/vpattern/complete/pattern4.moka"), mark_4template_volume);
  
    // TODO: Chargé séparément pour le moment, voir si je laisse comme ça 
    partial_3_template.load_additional_vpattern(CGAL::data_file_path("hexmeshing/vpattern/partial/pattern3.moka"), mark_3template_partial_volume);
  }

  void load_surface(const std::string& file, Polyhedron& out_surface, Tree& out_aab) {
    std::ifstream off_file(file);
    CGAL_precondition_msg(off_file.good(), ("Input .off couldn't be read : " + file).c_str());

    Polyhedron surface;
    off_file>>out_surface;
    CGAL::Polygon_mesh_processing::triangulate_faces(out_surface);

    // Compute AABB tree
    out_aab.insert(faces(out_surface).first, faces(out_surface).second, out_surface);
    out_aab.accelerate_distance_queries();
    out_aab.bbox();
  }

  bool is_intersect(double x1, double y1, double z1,
                  double x2, double y2, double z2,
                  double x3, double y3, double z3,
                  double x4, double y4, double z4,
                  const Tree& t)
  {
    Kernel::Point_3 p1(x1,y1,z1);
    Kernel::Point_3 p2(x2,y2,z2);
    Kernel::Point_3 p3(x3,y3,z3);
    Kernel::Point_3 p4(x4,y4,z4);

    // And compute the two triangles
    Triangle t1(p1, p2, p3);
    if(t.do_intersect(t1))
    { return true; }

    t1=Triangle(p1, p3, p4);
    if(t.do_intersect(t1))
    { return true; }

    return false;
  }

  bool is_intersect(double x1, double y1, double z1,
                  double x2, double y2, double z2,
                  const Tree& t)
  {
    return
        is_intersect(x1,y1,z1, x2,y1,z1, x2,y1,z2, x1,y1,z2, t) || // f1 y1
        is_intersect(x2,y2,z1, x2,y1,z1, x2,y1,z2, x2,y2,z2, t) || // f2 x2
        is_intersect(x1,y2,z1, x2,y2,z1, x2,y2,z2, x1,y2,z2, t) || // f3 y2
        is_intersect(x1,y1,z1, x1,y1,z2, x1,y2,z2, x1,y2,z1, t) || // f4 x1
        is_intersect(x1,y1,z1, x2,y1,z1, x2,y2,z1, x1,y2,z1, t) || // f5 z1
        is_intersect(x1,y1,z2, x2,y1,z2, x2,y2,z2, x1,y2,z2, t);   // f6 z2
  }

  bool is_intersect(LCC& lcc, Dart_handle dh, const Tree& t)
  {
    CGAL::Bbox_3 bbox=lcc.point(dh).bbox();
    // For each vertex of the volume
    for(auto it=lcc.one_dart_per_incident_cell<0,3>(dh).begin(),
        itend=lcc.one_dart_per_incident_cell<0,3>(dh).end(); it!=itend; ++it)
    { bbox+=lcc.point(it).bbox(); }

    return is_intersect(bbox.xmin(), bbox.ymin(), bbox.zmin(),
                        bbox.xmax(), bbox.ymax(), bbox.zmax(), t);
  }

  void mark_intersecting_volume_with_poly(LCC& lcc, Polyhedron& poly, Tree& tree, Dart_handle dart) {
    if (is_intersect(lcc, dart, tree)) {
      get_or_create_attr<3>(lcc, dart);
    }
  }
} 

namespace CGAL::HexRefinement {

  // Identifies which 3-cell should be refined 
  using MarkingFunction = std::function<void(LCC&, Polyhedron&, Tree&, Dart_handle)>;

  LCC tworefinement(const std::string& file, int nb_subdiv_per_dim, MarkingFunction cellIdentifier) {
    using namespace TwoRefinement;

    HexMeshingData hdata;
    LCC& lcc = hdata.lcc;
    load_patterns(hdata.regular_templates, hdata.partial_templates);

    Polyhedron surface; 
    Tree aabb;
    load_surface(file, surface, aabb);
  
    hdata.grid = generate_grid(hdata.lcc, aabb, nb_subdiv_per_dim);
    hdata.identified_mark = lcc.get_new_mark();
    hdata.template_mark = lcc.get_new_mark();
    hdata.corner_mark = lcc.get_new_mark();

    lcc.negate_mark(hdata.corner_mark);
    
    for (auto dart = lcc.one_dart_per_cell<3>().begin(), end = lcc.one_dart_per_cell<3>().end(); dart != end; dart++){
      cellIdentifier(lcc, surface, aabb, dart);
    }

    extract_first_faces_of_planes(hdata);

    for (int p = 0; p < 1; p++) {
      PlaneData pdata;

      mark_identified_cells_from_3_attrs(hdata);

      extract_darts_from_even_planes(hdata, pdata, (Plane)p);

      if (p > 0){
        propagation_stage(hdata, pdata);
      }

      // fix_impossible_cases(hdata, pdata);

      create_vertices_for_templates(hdata, pdata);

      refine_marked_hexes(hdata, pdata);
    }

    debug_node_mark = hdata.template_mark;
    // lcc.free_mark(corner_mark);
    return lcc;
  } 

}
