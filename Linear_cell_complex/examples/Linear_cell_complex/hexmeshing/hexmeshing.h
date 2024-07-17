#pragma once

#include <CGAL/Graphics_scene.h>
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
#include <bitset>
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
  struct FaceAttrValue { char template_id = 0; std::bitset<3> plane; std::bitset<3> explored; };
  struct VolumeAttrValue { char iteration = -1; bool draw = false; };

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
size_type debug3;
std::vector<Dart_handle> suspect;
namespace CGAL::HexRefinement::TwoRefinement {
  enum Plane { YZ, ZX, XY, YX = XY, ZY = YZ, XZ = ZX };


  struct Grid {
    Point center;
    double size;
    int nb_subdiv_per_dim;
    Grid(): center(0,0,0), size(0), nb_subdiv_per_dim(0){}
    Grid(Point c, double s, int n): center(c), size(s), nb_subdiv_per_dim(n){}
  };

  struct HexMeshingData {
    LCC lcc;
    size_type identified_mark, template_mark, corner_mark, propagation_face_mark;
    Pattern_substituer<LCC> regular_templates, partial_templates;
    Grid grid;

    std::array<std::vector<Dart_handle>, 3> first_face_of_planes;
  };

  struct PlaneData {
    Plane plane;

    std::vector<Dart_handle> marked_nodes;
    std::vector<Dart_handle> additionnal_volumes_found;

    std::vector<Dart_handle> volumes_to_refine;
    std::vector<Dart_handle> faces_to_refine; // Array =  [ faces of even planes , faces adjacent to a marked edege of even planes ]
    size_t face_of_even_planes_end = -1;

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

  void mark_face_unchecked(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      lcc.mark(dit, mark);

      if (!lcc.is_free<3>(dit))
        lcc.mark(lcc.beta<3>(dit), mark);
    }

    assert(lcc.is_whole_cell_marked<2>(dart, mark));
  }

  bool is_half_face_marked(LCC& lcc, Dart_handle dart, size_type mark ){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      if (!lcc.is_marked(dit, mark)) return false;
    }

    return true;
  }

  void mark_half_face_unchecked(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);

    // // DEBUG
    // bool marked_other_side = false;
    // int c = 0;
    // for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
    //   if (!lcc.is_free<3>(dit) && lcc.is_marked(lcc.beta(dit, 3), mark)) {
    //     c++;
    //   }
    // }

    // assert(c == 4 || c == 0);

    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      lcc.mark(dit, mark);
    }

    assert(is_half_face_marked(lcc, dart, mark));
    // if (c == 0 && !lcc.is_free<3>(dart)) assert(!is_half_face_marked(lcc, lcc.beta(dart, 3), mark));
  }

  // Not the right way to do it, i don't want to have two copy of the same function and just one line changes
  // If adjacent face doesn't exist, returns lcc.null_dart_descriptor.
  template <bool add_incident_volumes = false>
  Dart_handle adjacent_face_on_plane(LCC& lcc, PlaneData &pdata, Dart_handle edge){
    Dart_handle other_face = lcc.beta(edge, 2);
    auto other_face_handle = lcc.attribute<2>(other_face);

    // Early exit if the beta(2) face is already on the plane
    if (other_face_handle != nullptr && other_face_handle->info().plane[pdata.plane])
      return other_face;

    // Early exit if no volume exists next to the beta(2) face
    if (lcc.is_free<3>(other_face))
      return lcc.null_dart_descriptor;

    bool found = false;
    for (int i = 0; i < 5 && !found; i++){
      // there is at most 3 steps to get to the face normal to the plane
      // if we are still in the for loop after 3 steps => bug
      assert(i < 4);

      // there must be no additionnal volumes on the first iteration.
      assert( ! (pdata.plane == 0 && i > 0) );

      if (lcc.is_free<3>(other_face)) {
        // We hit the grid border, there are no adjacent face possible to that edge
        return lcc.null_dart_descriptor;
      }

      other_face = lcc.beta(other_face, 3, 2);
      other_face_handle = lcc.attribute<2>(other_face);
      found = other_face_handle != nullptr && other_face_handle->info().plane[pdata.plane];

      if (!found && add_incident_volumes){
        pdata.additionnal_volumes_found.push_back(other_face);
      }
    }

    assert(found);

    return other_face;
  }

  std::vector<Dart_handle> incident_faces_to_0_cell_on_plane(LCC& lcc, PlaneData &pdata, Dart_handle edge){
    std::vector<Dart_handle> arr;
    arr.reserve(4);

    arr.push_back(edge);

    Dart_handle d2 = lcc.beta(edge, 0);

    // Left neighbour
    Dart_handle adjacent_face = adjacent_face_on_plane(lcc, pdata, d2);
    bool left_neigh = adjacent_face != lcc.null_dart_descriptor;
    if (left_neigh) arr.push_back(adjacent_face);

    //Right neighbour
    adjacent_face = adjacent_face_on_plane(lcc, pdata, edge);
    bool right_neigh = adjacent_face != lcc.null_dart_descriptor;
    if (right_neigh) arr.push_back(lcc.beta(adjacent_face, 1));

    // Diagonal neighbour to 'edge' face
    if (left_neigh && right_neigh){
      Dart_handle d3 = adjacent_face_on_plane(lcc, pdata, lcc.beta(arr[1], 0));
      if (d3 != lcc.null_dart_descriptor){
        arr.push_back(d3);
        assert(arr[2] == adjacent_face_on_plane(lcc, pdata, lcc.beta(d3, 0)));
      }
    }

    return arr;
  }

  void clean_up_3_template(HexMeshingData &hdata, const Dart_handle &origin_dart, const Dart_handle &upper_edge, const Dart_handle lower_edge, Dart_handle &face1, Dart_handle &face2)
  {
    LCC& lcc = hdata.lcc;

    // Ne pas maintenir des refs sur les attributs, ils vont disparaitre
    DartInfo::FaceAttrValue face1_attr, face2_attr;
    if (lcc.attribute<2>(face1) != nullptr) face1_attr = lcc.attribute<2>(face1)->info();
    if (lcc.attribute<2>(face2) != nullptr) face2_attr = lcc.attribute<2>(face2)->info();


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

    // Requires at least one face for merging/adding 2-attributes
    if (face1 == nullptr && face2 == nullptr ) return;

    // Merge the two previous face attributes bitsets
    std::bitset<3> merged_planes = face1_attr.plane | face2_attr.plane;

    auto merge_face = face1 != nullptr ? face1 : face2;
    auto merge_face_handle = get_or_create_attr<2>(lcc, merge_face);

    merge_face_handle->info().plane = merged_planes;
    suspect.push_back(merge_face);
  }

  template <unsigned int i, typename... DartArray>
  void assert_dart_attr_are_unique(LCC& lcc, DartArray... array){
    std::unordered_map<void*, int> attributes;

    ([&]{
      int index = 0;
      for (auto dart : array){
        auto attr = lcc.attribute<i>(dart);
        assert(attr != nullptr);
        auto ptr = &attr->info();
        assert(attributes.count(ptr) == 0);
        attributes[ptr] = index;
        index++;
      }
    }(),...);
  }

  void assert_faces_of_plane_valid(HexMeshingData& hdata, PlaneData& pdata){
    LCC& lcc = hdata.lcc;
    assert_dart_attr_are_unique<2>(hdata.lcc, pdata.faces_to_refine);
    
    for (int i = 0; i < pdata.face_of_even_planes_end; i++){
      Dart_handle face = pdata.faces_to_refine[i];
      auto& attr = lcc.attribute<2>(face)->info();
      auto nodes = lcc.darts_of_cell<2, 0>(face);
      int marked = 0;
      for (auto it = nodes.begin(), end = nodes.end(); it != end; it++){
        if (lcc.is_marked(it, hdata.template_mark)) marked++;
      }

      assert( attr.explored[pdata.plane] );
      assert( attr.plane[pdata.plane]  == true );
      assert( attr.template_id == marked );
    }
  }

  // void assert_templates_are_valid(HexMeshingData &hdata, PlaneData &pdata, std::vector<Dart_handle> array){
  //   LCC& lcc = hdata.lcc;
  //   int i = 0;
  //   for (auto vol : array){
  //     auto &attr = lcc.attribute<3>(vol)->info();

  //     int count = 0;
  //     auto d = lcc.darts_of_cell<2,0>(vol);
  //     for (auto it = d.begin(); it != d.end(); it++){
  //       if (lcc.is_marked(it, hdata.template_mark)) count++;
  //     }

  //     auto d2 = lcc.darts_of_cell<2,0>(lcc.beta(vol, 2, 1, 1, 2));
  //     for (auto it = d2.begin(); it != d2.end(); it++){
  //       if (lcc.is_marked(it, hdata.template_mark)) count++;
  //     }

  //     // assert( count == attr.template_id );
  //     i++;
  //     if (count != attr.template_id){
  //       lcc.unmark_all(debug_edge_mark);
  //       mark_all_0_cells<3>(lcc, vol, debug_edge_mark);
  //       render(lcc, hdata.template_mark, debug_edge_mark);

  //       assert(false);
  //     }
  //   }
  // }

  // Ensures indirectly that all possible volumes are refined, since we know for sure create_vertices is valid
  // This is to ensure that we missed no volumes / faces for refinement.
  void assert_all_faces_are_quadrilateral(LCC &lcc){
    auto iterator = lcc.one_dart_per_cell<2>();
    lcc.unmark_all(debug_edge_mark);
    for (auto it = iterator.begin(); it != iterator.end(); it++){
      auto edges = lcc.darts_of_cell<2, 0>(it);
      // assert(edges.size() == 4);
      if (edges.size() != 4)
        mark_face_unchecked(lcc, it, debug_edge_mark);
    }
  }

  void assert_all_volumes_are_hexes(LCC &lcc){
    auto iterator = lcc.one_dart_per_cell<3>();
    for (auto it = iterator.begin(); it != iterator.end(); it++){
      auto edges = lcc.darts_of_cell<3, 0>(it);
      // assert(edges.size() == (4 * 6));
    }
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

    for (auto& dart : pdata.faces_to_refine)
    {
      // Utile uniquement si les faces marqués ne sont pas 100% templatés
      // auto edges_range = lcc.darts_of_cell<2, 1>(dart);
      // int marked_propagation = 0;
      // std::vector<Dart_handle> edges_vec;
      // for (auto it = edges_range.begin(), end = edges_range.end(); it != end; it++){
      //   edges_vec.push_back(it);
      //   if (lcc.is_marked(it, hdata.propagation_face_mark)) marked_propagation++;
      // }

      if (hdata.regular_templates.query_replace_one_face(lcc, dart, hdata.template_mark) != SIZE_T_MAX) nbsub++;

      // if (marked_propagation >= 4){
      //   for (auto dart : edges_vec){
      //     mark_half_face_unchecked(lcc, dart, hdata.propagation_face_mark);
      //   }
      // }
    }

    // Cannot easily assert if all faces has been correctly treated, because some faces don't have attr
    // and we don't refine 3/4 template faces.

    std::cout << nbsub << " face substitution was made" << std::endl;

    nbsub = 0;
    for (auto& dart : pdata.volumes_to_refine)
    {
      size_type temp = hdata.regular_templates.query_replace_one_volume(lcc, dart, hdata.template_mark);
      if (temp != SIZE_T_MAX) nbsub++;

      if (temp == 2) temp = 3;
      if (temp < 10) temp++;

      // assert(vol_attr.template_id >= 0 && vol_attr.template_id <= 4 && vol_attr.template_id != 3 || vol_attr.template_id == 8);
      // assert(vol_attr.template_id == temp || temp == SIZE_T_MAX && (vol_attr.template_id == 0 || vol_attr.template_id == 8));
    }

    // assert(nbsub == pdata.volumes_to_refine.size());

    // Refine remaining 3 patterns
    for (auto marked_face : pdata.partial_templates_to_refine){
      assert(lcc.attribute<2>(marked_face) != nullptr);
      assert(lcc.attribute<2>(marked_face)->info().template_id == 3);
      refine_3_template(hdata, marked_face);
      nbsub += 2;
    }

    std::cout << nbsub << " volumic substitution was made" << std::endl;
  }

  void __expand_0_cell_marking(LCC &lcc, PlaneData &pdata, std::queue<Dart_handle>& faces_to_check, Dart_handle &edge) {
    auto faces = incident_faces_to_0_cell_on_plane(lcc, pdata, edge);

    for (Dart_handle face : faces){
      assert( lcc.attribute<2>(face) != nullptr);
      auto& face_attr =  lcc.attribute<2>(face)->info();

      // If the face didn't have any template before, it will have one, so add it in faces to refine
      if (face_attr.template_id == 0)
        pdata.faces_to_refine.push_back(face);

      face_attr.template_id++;
      assert(face_attr.template_id <= 4);

      if (face_attr.template_id == 2 || face_attr.template_id == 3)
        faces_to_check.push(face);
    }
  }

  bool fix_adjacent_3_templates(HexMeshingData &hdata, PlaneData &pdata, std::queue<Dart_handle>& faces_to_check, Dart_handle &face)
  {
    // Transform nearby 3 templates into two 4 templates
    // if they both share the unmarked node

    LCC& lcc = hdata.lcc;

    auto &face_attr = lcc.attribute<2>(face)->info();

    if (face_attr.template_id != 3) return false;

    auto edges = lcc.darts_of_cell<2, 1>(face);

    // find the unmarked edge;
    bool found = false;
    Dart_handle edges_unmarked[2];

    for (auto it = edges.begin(); it != edges.end() && !found; it++)
    {
      if (!lcc.is_marked(it, hdata.template_mark)){
        found = true;
        edges_unmarked[0] = it;
      }
    }

    assert(found);
    edges_unmarked[1] = lcc.beta(edges_unmarked[0], 0);

    bool is_impossible = false;
    for (int i = 0; i < 2; i++){
      Dart_handle edge = edges_unmarked[i];
      Dart_handle other_face = adjacent_face_on_plane(lcc, pdata, edge);

      if (other_face == lcc.null_dart_descriptor) continue;

      auto other_face_handle = lcc.attribute<2>(other_face);

      assert(other_face_handle != nullptr);

      if (other_face_handle->info().template_id != 3)
        continue;

      is_impossible = true;
      break;
    }

    if (is_impossible) {
      lcc.mark_cell<0>(edges_unmarked[0], hdata.template_mark);
      pdata.marked_nodes.push_back(edges_unmarked[0]);
      __expand_0_cell_marking(lcc, pdata, faces_to_check, edges_unmarked[0]);

      return true;
    }

    return false;

  }

  // If not all marks are chained consecutively, mark the whole face
  bool fix_mark_connectivity(HexMeshingData &hdata, PlaneData &pdata, std::queue<Dart_handle>& faces_to_check, Dart_handle face){
    LCC& lcc = hdata.lcc;

    auto edges = lcc.darts_of_cell<2, 1>(face);

    bool connected = true;
    Dart_handle edge1 = face;
    if (lcc.is_marked(edge1, hdata.template_mark)
      && !lcc.is_marked(lcc.beta(edge1, 1), hdata.template_mark)
      &&  lcc.is_marked(lcc.beta(edge1, 1, 1), hdata.template_mark))
    { connected = false;}

    Dart_handle edge2 = lcc.beta(face, 1);
    if (lcc.is_marked(edge2, hdata.template_mark)
      && !lcc.is_marked(lcc.beta(edge2, 1), hdata.template_mark)
      &&  lcc.is_marked(lcc.beta(edge2, 1, 1), hdata.template_mark))
    { connected = false;}

    if (connected) return false;

    // mark the face
    for (auto edge = edges.begin(); edge != edges.end(); edge++) {
      if (lcc.is_marked(edge, hdata.template_mark)) continue;

      lcc.mark_cell<0>(edge, hdata.template_mark);
      pdata.marked_nodes.push_back(edge);

      __expand_0_cell_marking(lcc, pdata, faces_to_check, edge);
    }

    auto &face_attr = lcc.attribute<2>(face)->info();
    face_attr.template_id = 4;

    return true;

  }

  void fix_impossible_cases(HexMeshingData &hdata, PlaneData &pdata){
    int fix_c_count = 0, fix_3_count = 0;

    // Condition for checking faces : 2 or 3 templates.
    std::queue<Dart_handle> faces_to_check;

    // First iteration

    // face_to_refine will grow, we only want to examine only faces that existed right now
    int faces_end = pdata.faces_to_refine.size();

    for (int i = 0; i < faces_end; i++){
      Dart_handle face = pdata.faces_to_refine[i];
      auto& face_attr = hdata.lcc.attribute<2>(face)->info();

      if (face_attr.template_id == 2 && fix_mark_connectivity(hdata, pdata, faces_to_check, face))
        fix_c_count++;

      if (face_attr.template_id == 3 && fix_adjacent_3_templates(hdata, pdata, faces_to_check, face))
        fix_3_count++;
    }

    // Repeat until there are no more faces to check
    while (!faces_to_check.empty()){
      Dart_handle front_face = faces_to_check.front();
      faces_to_check.pop();

      auto& face_attr = hdata.lcc.attribute<2>(front_face)->info();

      if (face_attr.template_id == 2 && fix_mark_connectivity(hdata, pdata, faces_to_check, front_face))
        fix_c_count++;

      if (face_attr.template_id == 3 && fix_adjacent_3_templates(hdata, pdata, faces_to_check, front_face))
        fix_3_count++;
    }

    std::cout << "Diagonal 2 templates fixed: " << fix_c_count << std::endl;
    std::cout << "Neighboring 3 templates repaired: " << fix_3_count << std::endl;
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
  void explore_face_of_plane(HexMeshingData& hdata, PlaneData& pdata, std::queue<Dart_handle>& queue, Dart_handle face, size_type explored_mark, size_type explored_edge) {
    LCC& lcc = hdata.lcc;

    auto edges = lcc.darts_of_cell<2,1>(face);

    auto& face_attr = get_or_create_attr<2>(lcc, face)->info();

    if (face_attr.explored[pdata.plane]) return;

    // The 3-attr might not exist if the cell is not identified.
    face_attr.explored[pdata.plane] = true;
    face_attr.template_id = 0;

    int nb = 0;
    // Add neighboring faces
    for (auto dit = edges.begin(), dend = edges.end(); dit != dend; dit++){
      bool explored = lcc.is_marked(dit, explored_mark);
      bool edge_explored = lcc.is_whole_cell_marked<1>(dit, explored_edge);
      bool identified = lcc.is_marked(dit, hdata.identified_mark);

      if (!explored){
        lcc.mark_cell<0>(dit, explored_mark);
      }

      if (!explored && identified){
        lcc.mark_cell<0>(dit, hdata.template_mark);
        pdata.marked_nodes.push_back(dit);
      }

      if (identified){
        // The 3attr is created only if one vertex is identiifed/marked
        face_attr.template_id++;
      }

      if (!edge_explored) {
        lcc.mark_cell<1>(dit, explored_edge);

        // Also add incident volumes while iterating
        Dart_handle other_face = adjacent_face_on_plane<true>(lcc, pdata, dit);
        // Also do the other way around
        if (!lcc.is_free<3>(dit))
          adjacent_face_on_plane<true>(lcc, pdata, lcc.beta(dit, 3));


        if (other_face != lcc.null_dart_descriptor)
          queue.push(other_face);
      }

      ++nb;
    }

    assert(nb == 4);

    if (face_attr.template_id > 0) {
      pdata.faces_to_refine.push_back(face);
    }
  }

  void mark_template_for_propagation(HexMeshingData &hdata, Dart_handle face, DartInfo::FaceAttrValue &face_attr)
  {
    LCC& lcc = hdata.lcc;
    Dart_handle marked_edge;
    auto nodes = lcc.darts_of_cell<2, 0>(face);

    bool found = false;

    if (face_attr.template_id == 1)
    {
      for (auto it = nodes.begin(); it != nodes.end() && !found; it++)
      {
        if (lcc.is_marked(it, hdata.template_mark))
        {
          found = true;
          marked_edge = it;
        }
      }
    }

    else if (face_attr.template_id == 2)
    {
      for (auto it = nodes.begin(); it != nodes.end() && !found; it++)
      {
        if (lcc.is_marked(it, hdata.template_mark) && lcc.is_marked(lcc.other_extremity(it), hdata.template_mark))
        {
          found = true;
          marked_edge = it;
        }
      }
    }

    assert(found);

    Dart_handle face1 = lcc.beta(marked_edge, 2, 1, 1, 2);
    Dart_handle face2 = lcc.beta(face1, 0, 2);
    Dart_handle face3 = lcc.beta(face1, 0, 0, 2);

    mark_half_face_unchecked(lcc, face1, hdata.propagation_face_mark);
    mark_half_face_unchecked(lcc, face3, hdata.propagation_face_mark);

    if (face_attr.template_id == 1)
    {
      mark_half_face_unchecked(lcc, face2, hdata.propagation_face_mark);
    }
  }

  void propagate_face(HexMeshingData &hdata, PlaneData &pdata, const Dart_handle &face, size_type explored_node_mark, size_type explored_face_mark)
  {
    LCC& lcc = hdata.lcc;
    Dart_handle back_face = lcc.beta(face, 2, 1, 1, 2);

    mark_all_0_cells<2>(lcc, back_face, hdata.template_mark);
    assert(lcc.attribute<3>(face) != nullptr);

    auto &vol_attr = get_or_create_attr<3>(lcc, face)->info();
    vol_attr.iteration = pdata.plane;

    // /!\ ALSO: Mark add the hex attached to the back_face to the refinement.
    // It is garanted that faces and volumes added will be unique in our array
    // Because that back_hex is only accessible within the template itself, and not
    // accessible from the plane.

    Dart_handle back_volume_face = lcc.beta(back_face, 3);
    assert(back_volume_face != nullptr);

    auto &back_vol_attr = get_or_create_attr<3>(lcc, back_volume_face)->info();

    if (back_vol_attr.iteration != pdata.plane){
    pdata.volumes_to_refine.push_back(back_volume_face);
      back_vol_attr.iteration = pdata.plane;
    }

    // Also add the neighboring faces to that 4-template face forrefinement
    if (!lcc.is_marked(back_volume_face, explored_face_mark)){
    pdata.faces_to_refine.push_back(back_volume_face);
      mark_face_unchecked(lcc, back_volume_face, explored_face_mark);
    }

    auto edges = lcc.darts_of_cell<2, 1>(back_volume_face);

    for (auto it = edges.begin(); it != edges.end(); it++)
    {
      auto top_face = lcc.beta(it, 2);

      if (!lcc.is_marked(it, explored_node_mark)){
        pdata.marked_nodes.push_back(it);
        lcc.mark_cell<0>(it, explored_node_mark);
      }

      if (!lcc.is_whole_cell_marked<2>(top_face, explored_face_mark)){
        pdata.faces_to_refine.push_back(top_face);
        mark_face_unchecked(lcc, top_face, explored_face_mark);
      }
    }
  }

  void propagation_stage(HexMeshingData &hdata, PlaneData &pdata, size_type explored_node_mark, size_type explored_face_mark)
  {
    LCC& lcc = hdata.lcc;

    bool skip_propagation = pdata.plane == 0;
    int propagated_count = 0;
    int marked_for_prop_count = 0;

    std::vector<Dart_handle> faces_to_mark;

    for (int i = 0; i < pdata.face_of_even_planes_end; i++){
      Dart_handle face = pdata.faces_to_refine[i];
      auto& face_attr = lcc.attribute<2>(face)->info();

      // TODO attention propagation dans les volumes sandwhichés (?)

      // Propagate the face if possible
      if (!skip_propagation && face_attr.template_id == 4){
        if (is_half_face_marked(lcc, face, hdata.propagation_face_mark)){
          propagate_face(hdata, pdata, face, explored_node_mark, explored_face_mark);
          propagated_count++;
        }

        Dart_handle other_face = lcc.beta(face, 3);
        if (!lcc.is_free<3>(face) && is_half_face_marked(lcc, other_face, hdata.propagation_face_mark)) {
          propagate_face(hdata, pdata, other_face, explored_node_mark, explored_face_mark);
          propagated_count++;
        }
      }
      // Else mark the face for propagation if possible
      else if (face_attr.template_id == 1 || face_attr.template_id == 2){
        faces_to_mark.push_back(face);
      }
    }

    if (pdata.plane >= 2) return;

    lcc.unmark_all(hdata.propagation_face_mark);

    for (auto face : faces_to_mark){
        auto& face_attr = lcc.attribute<2>(face)->info();
        mark_template_for_propagation(hdata, face, face_attr);
        marked_for_prop_count++;

        if (!lcc.is_free<3>(face)){
          mark_template_for_propagation(hdata, lcc.beta(face, 3), face_attr);
          marked_for_prop_count++;
        }
    }

    std::cout << "Number of faces propagated : " << propagated_count << std::endl;
    std::cout << "Number of faces marked for propagation : " << marked_for_prop_count << std::endl;
  }

  void get_cells_to_refine_from_plane(HexMeshingData &hdata, PlaneData &pdata, size_type explored_faces)
  {
    LCC& lcc = hdata.lcc;

    // Iterate over all faces of even planes
    for (int i = 0; i < pdata.face_of_even_planes_end; i++)
    {
      Dart_handle face = pdata.faces_to_refine[i];

      auto edges = lcc.darts_of_cell<2, 0>(face);

      for (auto edge = edges.begin(); edge != edges.end(); edge++)
      {
        // get incident faces to marked node to be refined
        // Incident faces normal to the plane
        // We also need to prevent adding twice the faces by marking them

        if (!lcc.is_marked(edge, hdata.template_mark) && !lcc.is_marked(lcc.other_extremity(edge), hdata.template_mark))
          continue;

        auto top_face_1 = lcc.beta(edge, 2);
        auto top_face_2 = lcc.beta(edge, 3, 2);

        if (!lcc.is_whole_cell_marked<2>(top_face_1, explored_faces))
        {
          pdata.faces_to_refine.push_back(top_face_1);
          mark_face_unchecked(lcc, top_face_1, explored_faces);
        }

        if (top_face_2 != lcc.null_dart_descriptor && !lcc.is_whole_cell_marked<2>(top_face_2, explored_faces))
        {
          pdata.faces_to_refine.push_back(top_face_2);
          mark_face_unchecked(lcc, top_face_2, explored_faces);
        }
      }
      // Also add the adjacent volumes if there is atleast one marked node
      // Because we are on odd/even layers, we can't accidently add twice a volume

      auto &face_attr = lcc.attribute<2>(face)->info();
      auto &vol_attr = get_or_create_attr<3>(lcc, face)->info();

      bool vol_processed = vol_attr.iteration == pdata.plane;

      // Create first volume attr and push back the volume
      if (!vol_processed){
        vol_attr.iteration = pdata.plane;

        if (face_attr.template_id != 3){
          pdata.volumes_to_refine.push_back(face);
        }
        else
          pdata.partial_templates_to_refine.push_back(face); // Both volumes are treated if it is 3 template
      }

      // Create second volume attr and push the face from the other volume
      if (!lcc.is_free<3>(face))
      {
        Dart_handle other_face = lcc.beta(face, 3);
        auto &other_vol_attr = get_or_create_attr<3>(lcc, other_face)->info();
        bool other_vol_processed = other_vol_attr.iteration == pdata.plane;

        if (!other_vol_processed)
        {
          other_vol_attr.iteration = pdata.plane;

          if (face_attr.template_id != 3){
            pdata.volumes_to_refine.push_back(other_face);
          }
        }
      }
    }
  }

  void get_cells_to_refine_from_additionnal_volumes(HexMeshingData &hdata, PlaneData &pdata, size_type explored_face)
  {
    LCC& lcc = hdata.lcc;

    for (Dart_handle initial_edge : pdata.additionnal_volumes_found)
    {
      auto &vol_attr = lcc.attribute<3>(initial_edge)->info();

      // TODO We can accidentally twice (or thrice) the same volume
      if (vol_attr.iteration == pdata.plane)
        continue;

      vol_attr.iteration = pdata.plane;

      bool node_1_marked = lcc.is_marked(initial_edge, hdata.template_mark);
      bool node_2_marked = lcc.is_marked(lcc.other_extremity(initial_edge), hdata.template_mark);

      if (!node_1_marked && !node_2_marked)
        continue;

      Dart_handle adjacent_faces[] = {
        initial_edge,
        lcc.beta(initial_edge, 2),

        lcc.beta(initial_edge, 0, 2),
        lcc.beta(initial_edge, 1, 2)
      };

      // Compliqué pour pas grand chose:
      // Pour éviter d'itérerer sur toutes les darts des faces,
      // On sait que seul les noeuds de initial_edge peuvent être marqués

      if (node_1_marked || node_2_marked)
      {
        pdata.volumes_to_refine.push_back(initial_edge);

        if (!lcc.is_whole_cell_marked<2>(adjacent_faces[0], explored_face))
        {
          mark_face_unchecked(lcc, adjacent_faces[0], explored_face);
          pdata.faces_to_refine.push_back(adjacent_faces[0]);
        }

        if (!lcc.is_whole_cell_marked<2>(adjacent_faces[1], explored_face))
        {
          mark_face_unchecked(lcc, adjacent_faces[1], explored_face);
          pdata.faces_to_refine.push_back(adjacent_faces[1]);
        }
      }

      if (node_1_marked && !lcc.is_whole_cell_marked<2>(adjacent_faces[2], explored_face))
      {
        mark_face_unchecked(lcc, adjacent_faces[2], explored_face);
        pdata.faces_to_refine.push_back(adjacent_faces[2]);
      }

      if (node_2_marked && !lcc.is_whole_cell_marked<2>(adjacent_faces[3], explored_face))
      {
        mark_face_unchecked(lcc, adjacent_faces[3], explored_face);
        pdata.faces_to_refine.push_back(adjacent_faces[3]);
      }
    }
  }


  void extract_darts_from_even_planes(HexMeshingData &hdata, PlaneData &pdata, Plane iterationPlane)
  {
    LCC& lcc = hdata.lcc;
    pdata.plane = iterationPlane;

    size_type explored_mark = lcc.get_new_mark();
    size_type explored_face = lcc.get_new_mark();

    for (auto start_plane : hdata.first_face_of_planes[iterationPlane]) {
      std::queue<Dart_handle> to_explore;
      to_explore.push(start_plane); // First face

      while (!to_explore.empty()) {
        Dart_handle front = to_explore.front();
        to_explore.pop();
        explore_face_of_plane(hdata, pdata, to_explore, front, explored_mark, explored_face);
      }
    }

    fix_impossible_cases(hdata, pdata);

    pdata.face_of_even_planes_end = pdata.faces_to_refine.size();

    assert_faces_of_plane_valid(hdata, pdata);

    propagation_stage(hdata, pdata, explored_mark, explored_face);

    get_cells_to_refine_from_plane(hdata, pdata, explored_face);

    // Treat the volumes that at least share an edge with the plane
    get_cells_to_refine_from_additionnal_volumes(hdata, pdata, explored_face);


    lcc.free_mark(explored_mark);
    lcc.free_mark(explored_face);
  }


  void explore_and_setup_faces(HexMeshingData& hdata, std::queue<Dart_handle> &queue, Plane iteration, Dart_handle face, size_type explored_edge, size_type explored_face, int debug){
    LCC& lcc = hdata.lcc;
    auto edges = lcc.darts_of_cell<2,1>(face);

    if (lcc.is_whole_cell_marked<2>(face, explored_face)) return;

    auto& face_attr = get_or_create_attr<2>(lcc, face)->info();

    // The 3-attr might not exist if the cell is not identified.
    int nb = 0;

    face_attr.plane[iteration] = true;

    // Add neighboring faces
    for (auto dit = edges.begin(), dend = edges.end(); dit != dend; dit++){
      bool edge_explored = lcc.is_whole_cell_marked<1>(dit, explored_edge);
      bool has_opposite_face = !lcc.is_free<3>(lcc.beta(dit, 2));

      if (has_opposite_face && !edge_explored) {
        Dart_handle other_face = lcc.beta(dit, 2, 3, 2);
        lcc.mark_cell<1>(dit, explored_edge);

        if (other_face != lcc.null_dart_descriptor)
          queue.push(other_face);
      }

      ++nb;
    }

    assert(nb == 4);

    mark_face_unchecked(lcc, face, explored_face);
  }

  // Gather the faces first, if we don't we will have issues later manually
  // Iterating to the next plane after a first refinement stage
  void setup_even_planes(HexMeshingData& hdata){
    LCC& lcc = hdata.lcc;
    int nb_subdiv = hdata.grid.nb_subdiv_per_dim;

    // Extract the first faces of each plane
    for (int p = 0; p < 3; p++){
      Plane plane = (Plane)p;

      auto& start_of_planes = hdata.first_face_of_planes[p];

      Dart_handle start_plane = __first_vertex_of_even_plane(hdata.lcc, plane);

      for (int z = 0; z < nb_subdiv / 2; z++){
        start_of_planes.push_back(lcc.beta(start_plane, 0, 2));
        start_plane = __next_even_plane(lcc, start_plane, plane);
      }
    }

    // Create attributes of all faces of all planes
    // To be able to iterate properly at later stages
    size_type explored_edge = lcc.get_new_mark();
    size_type explored_face = lcc.get_new_mark();

    for (int p = 0; p < 3; p++){
      Plane plane = (Plane)p;

      int debug = 0;
      for (auto start_plane : hdata.first_face_of_planes[p]) {
        std::queue<Dart_handle> to_explore;
        to_explore.push(start_plane); // First face

        while (!to_explore.empty()) {
          Dart_handle front = to_explore.front();
          to_explore.pop();
          explore_and_setup_faces(hdata, to_explore, plane, front, explored_edge, explored_face, debug);
        }

        debug++;
      }

      lcc.unmark_all(explored_edge);
    }

    lcc.free_mark(explored_edge);
    lcc.free_mark(explored_face);

  }

  void mark_identified_cells_from_3_attrs(HexMeshingData& hdata) {
    LCC& lcc = hdata.lcc;

    auto identified = lcc.attributes<3>();

    for (auto it = identified.begin(), end = identified.end(); it != end; it++){
      mark_all_0_cells<3>(lcc, it->dart(), hdata.identified_mark);
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
    return Grid(pos, size, nb_subdiv_per_dim);
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

    int vertices_created = 0;
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

        vertices_created++;
        edges_to_subdivide.push_back(nit);
        lcc.mark_cell<1>(nit, arrete_done);
      }

      dart = lcc.beta<1>(dart);
    }

    for (Dart_handle dart : edges_to_subdivide)
    {
      lcc.insert_barycenter_in_cell<1>(dart);
    }

    std::cout << "Vertices created: " << vertices_created << std::endl;

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

  void render_two_refinement_result(const LCC& lcc, const char* title = "TwoRefinement Result"){
    LCCSceneOptions<LCC> gso;

    gso.draw_volume = [](const LCC& lcc, LCC::Dart_const_handle dart){return lcc.attribute<3>(dart) != nullptr;};
    gso.colored_volume = [](const LCC& lcc, LCC::Dart_const_handle dart){ return true; };
    gso.volume_color = [](const LCC& lcc, LCC::Dart_const_handle dart){ return rand_color_from_dart(lcc, dart); };

    CGAL::Graphics_scene buffer;
    add_to_graphics_scene(lcc, buffer, gso);
    CGAL::draw_graphics_scene(buffer);
  }

  void debug_render(const LCC& lcc, TwoRefinement::HexMeshingData& hdata, TwoRefinement::PlaneData& pdata, const char* title = "TwoRefinement Result"){
    LCCSceneOptions<LCC> gso;

    gso.draw_volume = [](const LCC& lcc, LCC::Dart_const_handle dart){return lcc.attribute<3>(dart) != nullptr;};
    // gso.volume_color = [](const LCC& lcc, LCC::Dart_const_handle dart){ return red(); };
    // gso.colored_volume = [](const LCC& lcc, LCC::Dart_const_handle dart){ return lcc.attribute<3>(dart) != nullptr && lcc.attribute<3>(dart)->info().draw; };

    // gso.colored_edge = [](const LCC& lcc, LCC::Dart_const_handle dart){ return lcc.is_whole_cell_marked<1>(dart, debug3); };
    // gso.edge_color = [](const LCC& lcc, LCC::Dart_const_handle dart){ return green(); }

    // ============================= WORKS
    // gso.draw_face = [](const LCC& lcc, LCC::Dart_const_handle dart){return true;};
    // gso.face_color = [](const LCC& lcc, LCC::Dart_const_handle dart){
    //   auto a = lcc.attribute<2>(dart);
    //   return a != nullptr && a->info().plane[2] ? red() : lcc.is_whole_cell_marked<2>(dart, debug_edge_mark) ? green() : white();};
    // gso.colored_face = [&](const LCC& lcc, LCC::Dart_const_handle dart){return true;};

    // gso.colored_vertex = [](const LCC& lcc, LCC::Dart_const_handle dart){return true;};
    // gso.vertex_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){return lcc.is_marked(dart, hdata.template_mark) ? green() : blue();};
    // =============================

    gso.draw_volume = [](const LCC& lcc, LCC::Dart_const_handle dart){
      return true;
    };
    // gso.colored_volume = [](const LCC& lcc, LCC::Dart_const_handle dart){
      
    // };
    // gso.volume_color = [](const LCC& lcc, LCC::Dart_const_handle dart){
    //   return red();
    // };

    // gso.face_wireframe = [](const LCC& lcc, LCC::Dart_const_handle dart){return true;};

    gso.draw_face = [](const LCC& lcc, LCC::Dart_const_handle dart){return true;};
    gso.face_color = [](const LCC& lcc, LCC::Dart_const_handle dart){
      auto a = lcc.attribute<2>(dart);
      auto b = lcc.attribute<3>(dart);
      return a != nullptr && a->info().plane[2] ? red() : lcc.is_whole_cell_marked<2>(dart, debug_edge_mark) ? green() : white();};
    gso.colored_face = [&](const LCC& lcc, LCC::Dart_const_handle dart){return true;};

    gso.colored_vertex = [](const LCC& lcc, LCC::Dart_const_handle dart){return true;};
    gso.vertex_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){return lcc.is_marked(dart, hdata.template_mark) ? green() : blue();};


    gso.edge_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){return green(); };
    gso.colored_edge = [&](const LCC& lcc, LCC::Dart_const_handle dart){return lcc.is_marked(dart, hdata.propagation_face_mark); };


    CGAL::Graphics_scene buffer;
    add_to_graphics_scene(lcc, buffer, gso);
    CGAL::draw_graphics_scene(buffer);
  }

  LCC two_refinement(const std::string& file, int nb_subdiv_per_dim, MarkingFunction cellIdentifier) {
    using namespace TwoRefinement;

    HexMeshingData hdata;
    LCC& lcc = hdata.lcc;
    load_patterns(hdata.regular_templates, hdata.partial_templates);

    Polyhedron surface;
    Tree aabb;
    load_surface(file, surface, aabb);

    debug_node_mark = lcc.get_new_mark();
    debug_edge_mark = lcc.get_new_mark();

    hdata.grid = generate_grid(hdata.lcc, aabb, nb_subdiv_per_dim);
    hdata.identified_mark = lcc.get_new_mark();
    hdata.template_mark = lcc.get_new_mark();
    hdata.corner_mark = lcc.get_new_mark(); // TODO not useful to the algorithm, only for checks in 3template
    hdata.propagation_face_mark = lcc.get_new_mark();

    for (auto dart = lcc.one_dart_per_cell<3>().begin(), end = lcc.one_dart_per_cell<3>().end(); dart != end; dart++){
      cellIdentifier(lcc, surface, aabb, dart);
    }

    setup_even_planes(hdata);

    for (int p = 0; p < 3; p++) {
      PlaneData pdata;

      lcc.negate_mark(hdata.corner_mark);

      mark_identified_cells_from_3_attrs(hdata);

      extract_darts_from_even_planes(hdata, pdata, (Plane)p);

      assert_dart_attr_are_unique<3>(lcc, pdata.volumes_to_refine, pdata.partial_templates_to_refine);

      create_vertices_for_templates(hdata, pdata);

      refine_marked_hexes(hdata, pdata);

      // debug3 = lcc.get_new_mark();
      for (auto a : pdata.additionnal_volumes_found){

        if (pdata.plane != 2) continue;
        
        get_or_create_attr<3>(lcc, a)->info().draw = true;
      }

      for (auto a : pdata.faces_to_refine){

        if (pdata.plane != 2) continue;
        
        get_or_create_attr<3>(lcc, a)->info().draw = true;
      }

      assert_all_faces_are_quadrilateral(lcc);
      debug_render(lcc, hdata, pdata);
      // lcc.free_mark(debug3);

      assert_all_volumes_are_hexes(lcc);

      lcc.unmark_all(hdata.identified_mark);
      lcc.unmark_all(hdata.template_mark);
      lcc.unmark_all(hdata.corner_mark);


    }

    // debug_node_mark = hdata.template_mark;
    // lcc.free_mark(corner_mark);
    return lcc;
  }

}
