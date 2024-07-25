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
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/Union_find.h>
#include <CGAL/aff_transformation_tags.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <bitset>
#include <functional>
#include <guiddef.h>
#include "utils.h"

#include "../query_replace/cmap_query_replace.h"
#include "../query_replace/init_to_preserve_for_query_replace.h"

#include <boost/container/static_vector.hpp>
#include <limits>
#include <qregularexpression.h>
#include <type_traits>
#include <unordered_map>
#include <vector>

#define NDEBUG

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT       FT;
typedef Kernel::Point_3  Point;
typedef Kernel::Vector_3 Vector;

typedef typename Kernel::Triangle_3 Triangle;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_Traits;
typedef CGAL::AABB_tree<AABB_Traits> Tree;


enum class VolumeType {
  NONE,         // Newly created vol_attribute
  REFINEMENT,   // Previously NONE volumes that was selected for refinement
  ID_EXPANSION, // Expansion of the Identified layer of volumes
  IDENTIFIED    // Volumes identified for refinement
};

class LCCItems
{
public:
  template < class Storage >
  struct Dart_wrapper
  {
    struct VolumeAttrValue {
      char iteration = -1;
      VolumeType type = VolumeType::NONE;
    };

    struct FaceAttrValue {
      char template_id = 0;
      std::bitset<3> plane;
      std::bitset<3> explored;
      typename CGAL::Union_find<typename Storage::Dart_handle>::handle cc_id;
      bool even = false;
    };

    typedef CGAL::Cell_attribute<Storage, FaceAttrValue> FaceAttr;
    typedef CGAL::Cell_attribute<Storage, VolumeAttrValue> VolumeAttr;
    typedef std::tuple<CGAL::Cell_attribute_with_point<Storage>, void, FaceAttr, VolumeAttr> Attributes;
  };
};


typedef CGAL::Linear_cell_complex_traits<3,Kernel> LCCTraits;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3, LCCTraits, LCCItems> LCC;
typedef typename LCC::Dart_handle Dart_handle;
typedef typename LCC::Vertex_attribute_handle Vertex_handle;
typedef typename LCC::size_type  size_type;

using DartInfo = LCCItems::Dart_wrapper<LCC::Storage>;

// Do not use this, only to satisfy lcc.is_isomorphic_to
bool operator==(const DartInfo::VolumeAttrValue& first, const DartInfo::VolumeAttrValue& second){
  return true;
}

// Do no use this, only to statisfy lcc.is_isomorphic_to
bool operator==(const DartInfo::FaceAttrValue& first, const DartInfo::FaceAttrValue& second){
  return true;
}
const size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

// to be removed
size_type debug_node_mark;
size_type debug_edge_mark;
size_type debug3;
namespace CGAL::HexRefinement::TwoRefinement {
  enum Plane { YZ, ZX, XY, YX = XY, ZY = YZ, XZ = ZX };

  // Identifies which 3-cell should be refined
  using MarkingFunction = std::function<bool(LCC&, Polyhedron&, Tree&, Dart_handle)>;

  struct Grid {
    Point center;
    double size;
    int nb_subdiv_per_dim;
    Grid(): center(0,0,0), size(0), nb_subdiv_per_dim(0){}
    Grid(Point c, double s, int n): center(c), size(s), nb_subdiv_per_dim(n){}
  };

  struct HexMeshingData {
    LCC lcc;
    size_type identified_mark, template_mark, propagation_face_mark;
    Pattern_substituer<LCC> regular_templates, partial_templates;
    Grid grid;

    int level = 0;

    using PlaneCC = std::vector<Dart_handle>; // One dart per face connected components
    using PlaneSet = std::vector<PlaneCC>; // A set of planes

    std::array<PlaneSet, 3> first_face_of_planes;
  };

  struct RefinementData {
    Plane iteration;

    std::vector<Dart_handle> marked_nodes;
    std::vector<Dart_handle> additionnal_volumes_found;

    std::vector<Dart_handle> volumes_to_refine;
    std::vector<Dart_handle> faces_to_refine; // Array =  [ faces of even planes , faces adjacent to a marked edege of even planes ]
    size_t face_of_even_planes_end = -1;

    std::vector<Dart_handle> partial_templates_to_refine;
  };

  template <unsigned int i>
  typename LCC::Attribute_descriptor<i>::type get_or_create_attr(LCC& lcc, Dart_handle dart){
    auto attr = lcc.attribute<i>(dart);

    if (attr == nullptr){
      attr = lcc.create_attribute<i>();
      lcc.set_attribute<i>(dart, attr);
    }
    return attr;
  }

  LCC::Attribute_descriptor<3>::type get_or_create_refinement_volume(LCC& lcc, Dart_handle dart){
    auto attr = lcc.attribute<3>(dart);

    if (attr == nullptr){
      attr = lcc.create_attribute<3>();

      lcc.set_attribute<3>(dart, attr);
    }

    // Previously NONE volumes are tagged as refined
    if (attr->info().type == VolumeType::NONE)
      attr->info().type = VolumeType::REFINEMENT;

    return attr;

  }

  boost::container::static_vector<Dart_handle, 27> cells_26_connectivity(LCC& lcc, Dart_handle dart, bool include_self_vol = false) {
    boost::container::static_vector<Dart_handle, 27> array;
    Dart_handle layers[3] = {
      dart,
      lcc.beta(dart, 3),
      lcc.beta(dart, 2, 1, 1, 2, 3)
    };

    for (int i = 0; i < 3; i++){
      Dart_handle mid_face_dart = layers[i];

      if (mid_face_dart == lcc.null_dart_descriptor) continue;

      if (include_self_vol or i != 0) array.push_back(mid_face_dart);

      auto edges = lcc.darts_of_cell<2, 1>(mid_face_dart);
      for (auto edge_it = edges.begin(), end = edges.end(); edge_it != end; edge_it++){

        if (!lcc.is_free<3>(lcc.beta(edge_it, 2))){
          Dart_handle other_face = lcc.beta(edge_it, 2, 3, 2);
          if (other_face != lcc.null_dart_descriptor)
            array.push_back(other_face);

          // If nesting ..
          if (!lcc.is_free<3>(lcc.beta(other_face, 1, 2))){
            Dart_handle side_face = lcc.beta(other_face, 1, 2, 3, 2);
            if (side_face != lcc.null_dart_descriptor)
              array.push_back(side_face);
          }
        }

      }
    }

    return array;
  }

  Dart_handle find_3_template_origin(LCC& lcc, Dart_handle marked_face, size_type template_mark) {

    auto _edges_count = lcc.darts_of_cell<2,1>(marked_face).size();
    assert(_edges_count == 6);

    Dart_handle dart = marked_face;

    // Get the origin dart : Find the two unmarked node on the face
    // since the 3 template is created by adjacent two 2-templates
    bool found = false;
    for (int i = 0; i < 6; i++)
    {
      if (!lcc.is_marked(dart, template_mark)
        && !lcc.is_marked(lcc.beta(dart, 1), template_mark)
        && !lcc.is_marked(lcc.beta(dart, 1, 1), template_mark))
      {
        assert(lcc.is_marked(lcc.beta(dart, 1, 1, 1), template_mark));
        assert(lcc.is_marked(lcc.beta(dart, 1, 1, 1, 1), template_mark));
        assert(lcc.is_marked(lcc.beta(dart, 1, 1, 1, 1, 1), template_mark));
        assert(lcc.beta(dart, 1, 1, 1, 1, 1, 1) == dart);
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

  bool is_half_face_marked(const LCC& lcc, LCC::Dart_const_handle dart, size_type mark ){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);
    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      if (!lcc.is_marked(dit, mark)) return false;
    }

    return true;
  }

  void mark_half_face_unchecked(LCC& lcc, Dart_handle dart, size_type mark){
    auto iterator = lcc.darts_of_cell<2, 1>(dart);

    for (auto dit = iterator.begin(), dend = iterator.end(); dit != dend; dit++) {
      lcc.mark(dit, mark);
    }

    assert(is_half_face_marked(lcc, dart, mark));
  }

  Dart_handle __adjacent_face_on_plane(LCC& lcc, Plane plane, Dart_handle edge){
    Dart_handle other_face = lcc.beta(edge, 2);
    auto this_face_handle = lcc.attribute<2>(edge);
    auto other_face_handle = lcc.attribute<2>(other_face);

    // Early exit if the beta(2) face is already on the plane
    if (other_face_handle != nullptr
      && other_face_handle->info().plane[plane])
      return other_face;

    bool found = false;
    for (int i = 0; i < 5 && !found; i++){

      if (lcc.is_free<3>(other_face)) {
        // We hit the grid border, there are no adjacent face possible to that edge
        return lcc.null_dart_descriptor;
      }

      other_face = lcc.beta(other_face, 3, 2);
      other_face_handle = lcc.attribute<2>(other_face);

      found = other_face_handle != nullptr
        && other_face_handle->info().plane[plane]
        && this_face_handle != other_face_handle;
    }

    return found ? other_face : lcc.null_dart_descriptor;
  }

  Dart_handle __adjacent_face_on_plane(LCC& lcc, Plane plane, Dart_handle edge, std::vector<Dart_handle>& additional_volumes){
    Dart_handle other_face = lcc.beta(edge, 2);
    auto this_face_handle = lcc.attribute<2>(edge);
    auto other_face_handle = lcc.attribute<2>(other_face);
    boost::container::static_vector<Dart_handle, 5> __additional_volumes;

    // Early exit if the beta(2) face is already on the plane
    if (other_face_handle != nullptr
      && other_face_handle->info().plane[plane])
      return other_face;

    bool found = false;
    for (int i = 0; i < 5 && !found; i++){

      if (lcc.is_free<3>(other_face)) {
        // We hit the grid border, there are no adjacent face possible to that edge
        return lcc.null_dart_descriptor;
      }

      other_face = lcc.beta(other_face, 3, 2);
      other_face_handle = lcc.attribute<2>(other_face);

      found = other_face_handle != nullptr
        && other_face_handle->info().plane[plane]
        && this_face_handle != other_face_handle;

      if (!found){
        __additional_volumes.push_back(other_face);
      }
    }

    if (found)
      additional_volumes.insert(additional_volumes.end(), __additional_volumes.begin(), __additional_volumes.end());

    return found ? other_face : lcc.null_dart_descriptor;
  }


  Dart_handle adjacent_face_on_plane(LCC& lcc, Plane plane, Dart_handle edge){
    Dart_handle other_face = __adjacent_face_on_plane(lcc, plane, edge);

    if (other_face == lcc.null_dart_descriptor && !lcc.is_free<3>(edge))
      other_face = __adjacent_face_on_plane(lcc, plane, lcc.beta<3>(edge));

    return other_face;
  }

  boost::container::static_vector<Dart_handle, 4> incident_faces_to_0_cell_on_plane(LCC& lcc, RefinementData &rdata, Dart_handle edge){
    boost::container::static_vector<Dart_handle, 4> arr;

    // Add initial face
    arr.push_back(edge);
    Dart_handle d2 = lcc.beta(edge, 0);

    // Left neighbour
    Dart_handle adjacent_face = adjacent_face_on_plane(lcc, rdata.iteration, d2);
    bool left_neigh = adjacent_face != lcc.null_dart_descriptor;
    if (left_neigh) arr.push_back(adjacent_face);

    //Right neighbour
    adjacent_face = adjacent_face_on_plane(lcc, rdata.iteration, edge);
    bool right_neigh = adjacent_face != lcc.null_dart_descriptor;
    if (right_neigh) arr.push_back(lcc.beta(adjacent_face, 1));

    // Diagonal neighbour to 'edge' face
    if (left_neigh && right_neigh){
      int f = lcc.belong_to_same_cell<0>(arr[1], arr[0]) ? 0 : 1;
      Dart_handle d3 = adjacent_face_on_plane(lcc, rdata.iteration, lcc.beta(arr[1], f));
      assert(d3 != lcc.null_dart_descriptor);
      arr.push_back(d3);
      // TODO assertion only valid if call to __adjacent_face, and for others aswell
      //assert(arr[2] == adjacent_face_on_plane(lcc, rdata.iteration, lcc.beta(d3, 0)));
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
    std::bitset<3> merged_planes = face1_attr.plane | face2_attr.plane;
    auto merge_face = face1 != lcc.null_dart_descriptor ? face1 : face2;

    if (merge_face == lcc.null_dart_descriptor) return;

    // Merge the two previous face attributes bitsets
    auto& merge_face_attr = get_or_create_attr<2>(lcc, merge_face)->info();
    merge_face_attr = DartInfo::FaceAttrValue();
    merge_face_attr.plane = merged_planes;
  }

  template <unsigned int i, typename... DartArray>
  void assert_dart_attr_are_unique(LCC& lcc, DartArray... array){
    #ifndef NDEBUG
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
    #endif
  }

  void assert_faces_of_plane_valid(HexMeshingData& hdata, RefinementData& rdata){
    #ifndef NDEBUG
      LCC& lcc = hdata.lcc;
      assert_dart_attr_are_unique<2>(hdata.lcc, rdata.faces_to_refine);

      for (int i = 0; i < rdata.face_of_even_planes_end; i++){
        Dart_handle face = rdata.faces_to_refine[i];
        auto& attr = lcc.attribute<2>(face)->info();
        auto nodes = lcc.darts_of_cell<2, 0>(face);
        int marked = 0;
        for (auto it = nodes.begin(), end = nodes.end(); it != end; it++){
          if (lcc.is_marked(it, hdata.template_mark)) marked++;
        }

        assert( attr.explored[rdata.iteration] );
        assert( attr.plane[rdata.iteration]  == true );
        assert( attr.template_id == marked );
      }
    #endif
  }

  // Ensures indirectly that all possible volumes are refined, since we know for sure create_vertices is valid
  // This is to ensure that we missed no volumes / faces for refinement.
  void assert_all_faces_are_quadrilateral(LCC &lcc){
    #ifndef NDEBUG
      auto iterator = lcc.one_dart_per_cell<2>();
      for (auto it = iterator.begin(); it != iterator.end(); it++){
        auto edges = lcc.darts_of_cell<2, 0>(it);
        assert(edges.size() == 4);
      }
    #endif
  }

  void assert_all_volumes_are_hexes(LCC &lcc){
    #ifndef NDEBUG
      auto iterator = lcc.one_dart_per_cell<3>();
      for (auto it = iterator.begin(); it != iterator.end(); it++){
        auto edges = lcc.darts_of_cell<3, 0>(it);
        assert(edges.size() == (4 * 6));
      }
    #endif
  }

  void refine_3_template(HexMeshingData &hdata, Dart_handle marked_face)
  {
    // TODO Might be written better
    LCC& lcc = hdata.lcc;

    Dart_handle origin_dart = find_3_template_origin(lcc, marked_face, hdata.template_mark);
    Dart_handle vol2_origin_dart = lcc.beta(origin_dart, 3);

    Dart_handle upper_d1 = origin_dart;
    Dart_handle upper_d2 = lcc.beta(origin_dart, 1, 1);

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

    // contract the shared edge between the two volumes
    lcc.contract_cell<1>(upper_mid_1);
    lcc.contract_cell<1>(upper_mid_2);

    clean_up_3_template(hdata, origin_dart, upper_edge, lower_edge, face1, face2);
    clean_up_3_template(hdata, vol2_origin_dart, vol2_upper_edge, vol2_lower_edge, vol2_face1, vol2_face2);
  }

  void refine_marked_hexes(HexMeshingData& hdata, RefinementData& rdata)
  {
    LCC& lcc = hdata.lcc;
    int nbsub = 0;
    int nb_3_tp = 0;


    for (auto& dart : rdata.faces_to_refine)
    {
      //Utile uniquement si les faces marqués ne sont pas 100% templatés
      auto edges_range = lcc.darts_of_cell<2, 1>(dart);
      int propagation = 0, beta3_propagation = 0;
      bool has_beta3 = !lcc.is_free<3>(dart);

      std::vector<Dart_handle> edges_vec;
      for (auto it = edges_range.begin(), end = edges_range.end(); it != end; it++){
        edges_vec.push_back(it);
        if (lcc.is_marked(it, hdata.propagation_face_mark)) propagation++;
        if (lcc.is_marked(lcc.beta<3>(it), hdata.propagation_face_mark)) beta3_propagation++;
      }

      if (hdata.regular_templates.query_replace_one_face(lcc, dart, hdata.template_mark) != SIZE_T_MAX) nbsub++;

      if (propagation >= 1){
        for (auto dart : edges_vec){
          mark_half_face_unchecked(lcc, dart, hdata.propagation_face_mark);
        }
      }

      if (beta3_propagation >= 1){
        for (auto dart : edges_vec){
          mark_half_face_unchecked(lcc, lcc.beta<3>(dart), hdata.propagation_face_mark);
        }
      }
    }

    // Cannot easily assert if all faces has been correctly treated, because some faces don't have attr
    // and we don't refine 3/4 template faces.

    std::cout << nbsub << " face substitution was made" << std::endl;

    nbsub = 0;
    for (auto& dart : rdata.volumes_to_refine)
    {
      size_type temp = hdata.regular_templates.query_replace_one_volume(lcc, dart, hdata.template_mark);
      if (temp != SIZE_T_MAX) nbsub++;

      if (temp == 2) temp = 3;
      if (temp < 10) temp++;

      // assert(vol_attr.template_id >= 0 && vol_attr.template_id <= 4 && vol_attr.template_id != 3 || vol_attr.template_id == 8);
      // assert(vol_attr.template_id == temp || temp == SIZE_T_MAX && (vol_attr.template_id == 0 || vol_attr.template_id == 8));
    }

    // assert(nbsub == rdata.volumes_to_refine.size());

    // Refine remaining 3 patterns
    for (auto marked_face : rdata.partial_templates_to_refine){
      assert(lcc.attribute<2>(marked_face) != nullptr);
      assert(lcc.attribute<2>(marked_face)->info().template_id == 3);
      refine_3_template(hdata, marked_face);
      nbsub += 2;
    }

    std::cout << nbsub << " volumic substitution was made" << std::endl;
  }

  void __expand_0_cell_marking(LCC &lcc, RefinementData &rdata, std::queue<Dart_handle>& faces_to_check, Dart_handle &edge) {
    auto faces = incident_faces_to_0_cell_on_plane(lcc, rdata, edge);

    for (Dart_handle face : faces){
      assert( lcc.attribute<2>(face) != nullptr);
      auto& face_attr =  lcc.attribute<2>(face)->info();

      // If the face didn't have any template before, it will have one, so add it in faces to refine
      if (face_attr.template_id == 0)
        rdata.faces_to_refine.push_back(face);

      face_attr.template_id++;
      assert(face_attr.template_id <= 4);

      if (face_attr.template_id == 2 || face_attr.template_id == 3)
        faces_to_check.push(face);
    }
  }

  bool fix_adjacent_3_templates(HexMeshingData &hdata, RefinementData &rdata, std::queue<Dart_handle>& faces_to_check, Dart_handle &face)
  {
    // Transform nearby 3 templates into two 4 templates
    // if they both share the unmarked node

    LCC& lcc = hdata.lcc;

    auto &face_attr = lcc.attribute<2>(face)->info();

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
      Dart_handle other_face = adjacent_face_on_plane(lcc, rdata.iteration, edge);

      if (other_face == lcc.null_dart_descriptor) continue;

      auto other_face_handle = lcc.attribute<2>(other_face);

      assert(other_face_handle != nullptr);

      if (other_face_handle->info().template_id == 3){
        is_impossible = true;
        break;
      }

    }

    if (is_impossible) {
      lcc.mark_cell<0>(edges_unmarked[0], hdata.template_mark);
      rdata.marked_nodes.push_back(edges_unmarked[0]);
      __expand_0_cell_marking(lcc, rdata, faces_to_check, edges_unmarked[0]);

      return true;
    }

    return false;
  }

  // If not all marks are chained consecutively, mark the whole face
  bool fix_mark_connectivity(HexMeshingData &hdata, RefinementData &rdata, std::queue<Dart_handle>& faces_to_check, Dart_handle face){
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
      rdata.marked_nodes.push_back(edge);
      __expand_0_cell_marking(lcc, rdata, faces_to_check, edge);
    }

    auto &face_attr = lcc.attribute<2>(face)->info();
    face_attr.template_id = 4;

    return true;

  }

  void fix_impossible_cases(HexMeshingData &hdata, RefinementData &rdata){
    int fix_c_count = 0, fix_3_count = 0;

    // Condition for checking faces : 2 or 3 templates.
    std::queue<Dart_handle> faces_to_check;

    // First iteration

    // face_to_refine will grow, we only want to examine only faces that existed right now
    int faces_end = rdata.faces_to_refine.size();

    for (int i = 0; i < faces_end; i++){
      Dart_handle face = rdata.faces_to_refine[i];
      auto& face_attr = hdata.lcc.attribute<2>(face)->info();

      if (face_attr.template_id == 2 && fix_mark_connectivity(hdata, rdata, faces_to_check, face))
        fix_c_count++;

      if (face_attr.template_id == 3 && fix_adjacent_3_templates(hdata, rdata, faces_to_check, face))
        fix_3_count++;
    }

    // Repeat until there are no more faces to check
    while (!faces_to_check.empty()){
      Dart_handle front_face = faces_to_check.front();
      faces_to_check.pop();

      auto& face_attr = hdata.lcc.attribute<2>(front_face)->info();

      if (face_attr.template_id == 2 && fix_mark_connectivity(hdata, rdata, faces_to_check, front_face))
        fix_c_count++;

      if (face_attr.template_id == 3 && fix_adjacent_3_templates(hdata, rdata, faces_to_check, front_face))
        fix_3_count++;
    }

    std::cout << "Diagonal 2 templates fixed: " << fix_c_count << std::endl;
    std::cout << "Neighboring 3 templates repaired: " << fix_3_count << std::endl;
  }


  /**
   * Mark 0-cells
   * Gather 2-cells adjacent to marked 0-cells
   * Gather 3-cells adjacents to marked 0-cells
   */
  void explore_face_of_plane(HexMeshingData& hdata, RefinementData& rdata, std::queue<Dart_handle>& queue, Dart_handle face, size_type explored_mark, size_type explored_edge) {
    LCC& lcc = hdata.lcc;

    auto& face_attr = lcc.attribute<2>(face)->info();

    if (face_attr.explored[rdata.iteration]) return;

    lcc.mark_cell<2>(face, debug3);

    face_attr.explored[rdata.iteration] = true;
    face_attr.template_id = 0;

    int nb = 0;
    auto edges = lcc.darts_of_cell<2,1>(face);
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
        rdata.marked_nodes.push_back(dit);
      }

      if (identified){
        // The 3attr is created only if one vertex is identiifed/marked
        face_attr.template_id++;
      }

      if (!edge_explored) {
        lcc.mark_cell<1>(dit, explored_edge);

        // Also add incident volumes while iterating
        Dart_handle other_face = __adjacent_face_on_plane(lcc, rdata.iteration, dit, rdata.additionnal_volumes_found);
        // Also do the other way around to catch volumes in the opposite dir
        if (!lcc.is_free<3>(dit)){
          Dart_handle other_face2 = __adjacent_face_on_plane(lcc, rdata.iteration, lcc.beta(dit, 3), rdata.additionnal_volumes_found);
          if (other_face == lcc.null_dart_descriptor) other_face  = other_face2;
        }

        if (other_face != lcc.null_dart_descriptor)
          queue.push(other_face);
      }

      ++nb;
    }

    assert(nb == 4);

    if (face_attr.template_id > 0) {
      rdata.faces_to_refine.push_back(face);
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

  void propagate_face(HexMeshingData &hdata, RefinementData &rdata, const Dart_handle &face, size_type explored_node_mark, size_type explored_face_mark)
  {
    LCC& lcc = hdata.lcc;
    Dart_handle back_face = lcc.beta(face, 2, 1, 1, 2);

    mark_all_0_cells<2>(lcc, back_face, hdata.template_mark);
    assert(lcc.attribute<3>(face) != nullptr);

    auto &vol_attr = get_or_create_refinement_volume(lcc, face)->info();
    vol_attr.iteration = rdata.iteration;
    // /!\ ALSO: Mark add the hex attached to the back_face for refinement.
    // It is garanted that faces and volumes added will be unique in our array
    // Because that back_hex is only accessible within the template itself, and not
    // accessible from the plane.

    Dart_handle back_volume_face = lcc.beta(back_face, 3);
    assert(back_volume_face != nullptr);

    auto &back_vol_attr = get_or_create_refinement_volume(lcc, back_volume_face)->info();

    if (back_vol_attr.iteration != rdata.iteration){
      rdata.volumes_to_refine.push_back(back_volume_face);
      back_vol_attr.iteration = rdata.iteration;
    }

    // Also add the neighboring faces to that 4-template face forrefinement
    if (!lcc.is_marked(back_volume_face, explored_face_mark)){
      rdata.faces_to_refine.push_back(back_volume_face);
      mark_face_unchecked(lcc, back_volume_face, explored_face_mark);
    }

    auto edges = lcc.darts_of_cell<2, 1>(back_volume_face);

    for (auto it = edges.begin(); it != edges.end(); it++)
    {
      auto top_face = lcc.beta(it, 2);

      if (!lcc.is_marked(it, explored_node_mark)){
        rdata.marked_nodes.push_back(it);
        lcc.mark_cell<0>(it, explored_node_mark);
      }

      if (!lcc.is_whole_cell_marked<2>(top_face, explored_face_mark)){
        rdata.faces_to_refine.push_back(top_face);
        mark_face_unchecked(lcc, top_face, explored_face_mark);
      }
    }
  }

  void propagation_stage(HexMeshingData &hdata, RefinementData &rdata, size_type explored_node_mark, size_type explored_face_mark)
  {
    LCC& lcc = hdata.lcc;

    bool skip_propagation = rdata.iteration == 0;
    int propagated_count = 0;
    int marked_for_prop_count = 0;

    std::vector<Dart_handle> faces_to_mark;

    for (int i = 0; i < rdata.face_of_even_planes_end; i++){
      Dart_handle face = rdata.faces_to_refine[i];
      auto& face_attr = lcc.attribute<2>(face)->info();

      // TODO attention propagation dans les volumes sandwhichés (?)

      // Propagate the face if possible
      if (!skip_propagation && face_attr.template_id == 4){
        if (is_half_face_marked(lcc, face, hdata.propagation_face_mark)){
          propagate_face(hdata, rdata, face, explored_node_mark, explored_face_mark);
          propagated_count++;
        }

        Dart_handle other_face = lcc.beta(face, 3);
        if (!lcc.is_free<3>(face) && is_half_face_marked(lcc, other_face, hdata.propagation_face_mark)) {
          propagate_face(hdata, rdata, other_face, explored_node_mark, explored_face_mark);
          propagated_count++;
        }
      }
      // Else mark the face for propagation if possible
      else if (face_attr.template_id == 1 || face_attr.template_id == 2){
        faces_to_mark.push_back(face);
      }
    }

    if (rdata.iteration >= 2) return;

    //lcc.unmark_all(hdata.propagation_face_mark);

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

  void get_cells_to_refine_from_plane(HexMeshingData &hdata, RefinementData &rdata, size_type explored_faces)
  {
    LCC& lcc = hdata.lcc;

    // Iterate over all faces of even planes
    for (int i = 0; i < rdata.face_of_even_planes_end; i++)
    {
      Dart_handle face = rdata.faces_to_refine[i];

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
          rdata.faces_to_refine.push_back(top_face_1);
          mark_face_unchecked(lcc, top_face_1, explored_faces);
        }

        if (top_face_2 != lcc.null_dart_descriptor && !lcc.is_whole_cell_marked<2>(top_face_2, explored_faces))
        {
          rdata.faces_to_refine.push_back(top_face_2);
          mark_face_unchecked(lcc, top_face_2, explored_faces);
        }
      }
      // Also add the adjacent volumes if there is atleast one marked node
      // Because we are on odd/even layers, we can't accidently add twice a volume

      auto &face_attr = lcc.attribute<2>(face)->info();
      auto &vol_attr = get_or_create_refinement_volume(lcc, face)->info();

      bool vol_processed = vol_attr.iteration == rdata.iteration;

      // Two faces point to the same volume. We are on a merged cell (caused by 3 templates)
      // if (vol_processed){
      //   mark_all_0_cells<3>(lcc, face, hdata.template_mark);
      // }

      // Create first volume attr and push back the volume
      if (!vol_processed) {
        vol_attr.iteration = rdata.iteration;

        if (face_attr.template_id != 3){
          rdata.volumes_to_refine.push_back(face);
        }
        else
          rdata.partial_templates_to_refine.push_back(face); // Both volumes are treated if it is 3 template
      }

      // Create second volume attr and push the face from the other volume
      if (!lcc.is_free<3>(face))
      {
        Dart_handle other_face = lcc.beta(face, 3);
        auto &other_vol_attr = get_or_create_refinement_volume(lcc, other_face)->info();
        bool other_vol_processed = other_vol_attr.iteration == rdata.iteration;

        // if (other_vol_processed)
          // mark_all_0_cells<3>(lcc, other_face, hdata.template_mark);
        if (!other_vol_processed) {
          other_vol_attr.iteration = rdata.iteration;

          if (face_attr.template_id != 3){
            rdata.volumes_to_refine.push_back(other_face);
          }
        }
      }
    }
  }

  void get_cells_to_refine_from_additionnal_volumes(HexMeshingData &hdata, RefinementData &rdata, size_type explored_face)
  {
    LCC& lcc = hdata.lcc;

    // No additionnal volumes should be found on the first iteration
    assert(rdata.iteration != 0 || rdata.iteration == 0 && rdata.additionnal_volumes_found.size() == 0);

    for (Dart_handle initial_edge : rdata.additionnal_volumes_found)
    {
      auto &vol_attr = lcc.attribute<3>(initial_edge)->info();

      // TODO We can accidentally twice (or thrice) the same volume
      if (vol_attr.iteration == rdata.iteration)
        continue;

      vol_attr.iteration = rdata.iteration;

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
        rdata.volumes_to_refine.push_back(initial_edge);

        if (!lcc.is_whole_cell_marked<2>(adjacent_faces[0], explored_face))
        {
          mark_face_unchecked(lcc, adjacent_faces[0], explored_face);
          rdata.faces_to_refine.push_back(adjacent_faces[0]);
        }

        if (!lcc.is_whole_cell_marked<2>(adjacent_faces[1], explored_face))
        {
          mark_face_unchecked(lcc, adjacent_faces[1], explored_face);
          rdata.faces_to_refine.push_back(adjacent_faces[1]);
        }
      }

      if (node_1_marked && !lcc.is_whole_cell_marked<2>(adjacent_faces[2], explored_face))
      {
        mark_face_unchecked(lcc, adjacent_faces[2], explored_face);
        rdata.faces_to_refine.push_back(adjacent_faces[2]);
      }

      if (node_2_marked && !lcc.is_whole_cell_marked<2>(adjacent_faces[3], explored_face))
      {
        mark_face_unchecked(lcc, adjacent_faces[3], explored_face);
        rdata.faces_to_refine.push_back(adjacent_faces[3]);
      }
    }
  }

  void extract_darts_from_even_planes(HexMeshingData &hdata, RefinementData &rdata, Plane iterationPlane)
  {
    LCC& lcc = hdata.lcc;
    rdata.iteration = iterationPlane;

    size_type explored_mark = lcc.get_new_mark();
    size_type explored_face = lcc.get_new_mark();

    HexMeshingData::PlaneSet& plane_set = hdata.first_face_of_planes[iterationPlane];

    // Explore all even planes
    lcc.unmark_all(debug3);
    for (int i = 1; i < plane_set.size(); i += 2) {
      std::queue<Dart_handle> to_explore;

      for (auto start : plane_set[i])
        to_explore.push(start);

      while (!to_explore.empty()) {
        Dart_handle front = to_explore.front();
        to_explore.pop();
        explore_face_of_plane(hdata, rdata, to_explore, front, explored_mark, explored_face);
      }
    }

    if (hdata.level == 1)
    {
      LCCSceneOptions<LCC> gso;

      gso.draw_face = [](const LCC& lcc, LCC::Dart_const_handle dart){
        return true;
      };
      gso.face_color = [](const LCC& lcc, LCC::Dart_const_handle dart){
        auto a = lcc.attribute<2>(dart);
        auto b = lcc.attribute<3>(dart);
        return lcc.is_whole_cell_marked<2>(dart, debug3) ? red()
        : b != nullptr && b->info().type > VolumeType::REFINEMENT ? green()
        : blue();

        // return lcc.is_whole_cell_marked<3>(dart, debug3)? red() : a != nullptr && a->info().type == VolumeType::IDENTIFIED ? green() : blue();
      };
      gso.colored_face = [](const LCC& lcc, LCC::Dart_const_handle dart){ return true; };

      // gso.colored_vertex = [](const LCC& lcc, LCC::Dart_const_handle dart){
      //   return true;
      // };
      // gso.draw_vertex = [](const LCC& lcc, LCC::Dart_const_handle dart){
      //   return true;
      // };
      // gso.vertex_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){
      //   return lcc.is_marked(dart, hdata.template_mark) ? red() : blue();
      // };
      CGAL::Graphics_scene buffer;
      add_to_graphics_scene(lcc, buffer, gso);
      CGAL::draw_graphics_scene(buffer);
    }

    rdata.face_of_even_planes_end = rdata.faces_to_refine.size();
    assert_faces_of_plane_valid(hdata, rdata);

    fix_impossible_cases(hdata, rdata);

    rdata.face_of_even_planes_end = rdata.faces_to_refine.size();
    assert_faces_of_plane_valid(hdata, rdata);

    propagation_stage(hdata, rdata, explored_mark, explored_face);

    get_cells_to_refine_from_plane(hdata, rdata, explored_face);
    get_cells_to_refine_from_additionnal_volumes(hdata, rdata, explored_face);

    lcc.free_mark(explored_mark);
    lcc.free_mark(explored_face);
  }


  void explore_and_setup_faces(HexMeshingData& hdata, std::queue<Dart_handle> &queue, Plane iteration, Dart_handle face, size_type explored_edge, size_type explored_face, bool even){
    LCC& lcc = hdata.lcc;
    auto edges = lcc.darts_of_cell<2,1>(face);

    if (lcc.is_whole_cell_marked<2>(face, explored_face)) return;

    auto& face_attr = get_or_create_attr<2>(lcc, face)->info();

    // The 3-attr might not exist if the cell is not identified.
    int nb = 0;

    face_attr.plane[iteration] = true;
    face_attr.even = even;

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

  void clean_up_and_reevaluate_attributes(HexMeshingData& hdata, MarkingFunction& cellIdentifier, Polyhedron& surface, Tree& aabb){
    // Erase all volume / faces attributes outside the domain,
    // Reset all volumes inside the domain, and reevaluate their identification status
    //
    // Since the refinement is uniform, we subdivided each original cell in 8 new cells.
    // This means we only need to reevaluate identified cells, if they are still identified or not.
    LCC& lcc = hdata.lcc;

    auto face_attributes = lcc.attributes<2>();
    auto vol_attributes = lcc.attributes<3>();

    std::vector<LCC::Attribute_descriptor<2>::type> faces_to_delete;
    std::vector<LCC::Attribute_descriptor<3>::type> volumes_to_delete;

    for (auto it = face_attributes.begin(), end = face_attributes.end(); it != end; it++){
      auto vol_handle = lcc.attribute<3>(it->dart());
      auto other_vol_handle = lcc.attribute<3>(lcc.beta(it->dart(), 3));

      if ((vol_handle == nullptr or vol_handle->info().type <= VolumeType::REFINEMENT)
        && (other_vol_handle == nullptr or other_vol_handle->info().type <= VolumeType::REFINEMENT))
        faces_to_delete.push_back(it);
    }

    for (auto it = vol_attributes.begin(), end = vol_attributes.end(); it != end; it++){
      DartInfo::VolumeAttrValue old_info = it->info();
      if (old_info.type <= VolumeType::REFINEMENT ){
        volumes_to_delete.push_back(it);
        continue;
      }

      it->info() = DartInfo::VolumeAttrValue();

      if (old_info.type == VolumeType::IDENTIFIED && cellIdentifier(lcc, surface, aabb, it->dart()))
        it->info().type = VolumeType::IDENTIFIED;
    }



    size_type size_face_b = lcc.attributes<2>().size();
    size_type size_vol_b = lcc.attributes<3>().size();

    for (auto attr_desc : volumes_to_delete){
      lcc.set_attribute<3>(attr_desc->dart(), nullptr);
      lcc.erase_attribute<3>(attr_desc);
    }

    for (auto attr_desc : faces_to_delete){
      lcc.set_attribute<2>(attr_desc->dart(), nullptr);
      lcc.erase_attribute<2>(attr_desc);
    }

    size_type size_face = lcc.attributes<2>().size();
    size_type size_vol = lcc.attributes<3>().size();
  }

  void expand_identified_cells(HexMeshingData& hdata, int current_lvl, int nb_levels){
    LCC& lcc = hdata.lcc;

    if (nb_levels == 0 or current_lvl >= nb_levels) return;

    // Calculate the totaling cells needed per level (height), with i ranging from 0 (lowest level) to n (highest level)
    // Very simple sequence : 4, 6, 6, 6, 6, ....               (height of n cells, cell sized from current i-subdivision)
    auto height_of_refinement_level = [](int i, int nb_levels) -> int{
      return i == 0 ? 4 : 6; // One dart per c (inside a planeonnected components of faces
    }; // A set of planes

    //wip
    for (int i = 0 ; i < 2; i++){
      auto vol_attributes = lcc.attributes<3>();
      for(auto it = vol_attributes.begin(), end = vol_attributes.end(); it != end; it++){
        if (it->info().type != VolumeType::IDENTIFIED) continue;

        for (auto vol : cells_26_connectivity(lcc, it->dart())){
          get_or_create_attr<3>(lcc, vol)->info().type = VolumeType::ID_EXPANSION;
        }
      }
    }

  }

  void setup_initial_planes(HexMeshingData& hdata){
    LCC& lcc = hdata.lcc;
    int nb_subdiv = hdata.grid.nb_subdiv_per_dim;

    auto __first_plane = [](LCC& lcc, Plane plane){
      switch (plane){
        case YZ: return lcc.beta(lcc.first_dart(), 0, 2);
        case XY: return lcc.first_dart();
        case ZX: return lcc.beta(lcc.first_dart(), 1, 2, 1);
      }

      CGAL_assertion_msg(false, "Unexpected value for plane");
      return lcc.null_dart_descriptor;
    };

    auto __next_plane = [](LCC& lcc, Dart_handle start_plane, Plane plane){
      return lcc.beta(start_plane, 1, 2, 3, 2, 1);
    };

    // Extract the first faces of each plane
    for (int p = 0; p < 3; p++){
      Plane plane = (Plane)p;

      auto& plane_set = hdata.first_face_of_planes[p];

      Dart_handle start_plane = __first_plane(hdata.lcc, plane);

      for (int z = 0; z < nb_subdiv; z++){
        std::vector<Dart_handle> starts;
        starts.push_back(lcc.beta(start_plane, 0, 2));

        start_plane = __next_plane(lcc, start_plane, plane);

        plane_set.push_back(std::move(starts));
      }
    }

    // Create attributes of all faces of all planes
    // To be able to iterate properly at later stages
    size_type explored_edge = lcc.get_new_mark();
    size_type explored_face = lcc.get_new_mark();

    for (int p = 0; p < 3; p++){
      Plane plane = (Plane)p;

      int debug = 0;

      int i = 1;
      for (auto& plane_cc :  hdata.first_face_of_planes[p]) {
        std::queue<Dart_handle> to_explore;

        for (auto start : plane_cc)
          to_explore.push(start);

        while (!to_explore.empty()) {
          Dart_handle front = to_explore.front();
          to_explore.pop();
          explore_and_setup_faces(hdata, to_explore, plane, front, explored_edge, explored_face, i % 2 == 0);
        }

        i++;
      }

      lcc.unmark_all(explored_edge);
    }

    lcc.free_mark(explored_face);
  }

  // Gather the faces first, if we don't we will have issues later manually
  // Iterating to the next plane after a first refinement stage
  void setup_initial(HexMeshingData& hdata, MarkingFunction& cellIdentifier, Polyhedron& surface, Tree& aabb){
    LCC& lcc = hdata.lcc;
    setup_initial_planes(hdata);

    // Mark initial identified cells
    for (auto dart = lcc.one_dart_per_cell<3>().begin(), end = lcc.one_dart_per_cell<3>().end(); dart != end; dart++){
      if (cellIdentifier(lcc, surface, aabb, dart))
        get_or_create_attr<3>(lcc, dart)->info().type = VolumeType::IDENTIFIED;
    }
  }

  void setup_next_level_face(HexMeshingData& hdata,
                            std::queue<Dart_handle>& to_explore,
                            Union_find<Dart_handle>& odd_union_find,
                            Union_find<Dart_handle>& even_union_find,
                            Dart_handle face,
                            Plane planeIteration,
                            int edge_mark, int face_mark){
    LCC& lcc = hdata.lcc;
    auto vol_handle = lcc.attribute<3>(face);
    auto beta3_vol_handle = lcc.attribute<3>(lcc.beta<3>(face));

    bool identified = vol_handle != nullptr && vol_handle->info().type > VolumeType::REFINEMENT;
    bool beta3_identified = beta3_vol_handle != nullptr && beta3_vol_handle->info().type > VolumeType::REFINEMENT;

    if (lcc.is_whole_cell_marked<2>(face, face_mark)) return;
    lcc.mark_cell<2>(face, face_mark);

    auto edges = lcc.darts_of_cell<2, 1>(face);
    Union_find<Dart_handle>::handle odd_cc_id, even_cc_id;

    Dart_handle back_face = lcc.beta(face, 2, 1, 1, 2);

    for (auto it = edges.begin(), end = edges.end(); it != end; it++){
      // Iterate
      assert(lcc.belong_to_same_cell<3>(it, face));

      Dart_handle adjacent_face = __adjacent_face_on_plane(lcc, planeIteration, it);

      if (adjacent_face == lcc.null_dart_descriptor) continue;

      if (lcc.is_whole_cell_unmarked<1>(it, edge_mark)){
        to_explore.push(adjacent_face);
        lcc.mark_cell<1>(it, edge_mark);
      }


      // Odd plane union find, only if a volume attached to the face is identified
      if (identified or beta3_identified){
        auto& adj_face_attr = lcc.attribute<2>(adjacent_face)->info();

        if (odd_cc_id != nullptr
          && adj_face_attr.cc_id != nullptr
          && !odd_union_find.same_set(odd_cc_id, adj_face_attr.cc_id)){
          odd_union_find.unify_sets(odd_cc_id, adj_face_attr.cc_id);
        }

        if (odd_cc_id == nullptr && adj_face_attr.cc_id != nullptr)
          odd_cc_id = odd_union_find.find(adj_face_attr.cc_id);
      }

      // Even plane union find, only if the studied volume is identified
      if (identified){
        auto back_face_handle = lcc.attribute<2>(lcc.beta(adjacent_face, 2, 1, 1, 2));
        // Skip back face with no attributes
        if (back_face_handle == nullptr) continue;
        auto& back_face_attr = back_face_handle->info();

        if (even_cc_id != nullptr
          && back_face_attr.cc_id != nullptr
          && !even_union_find.same_set(even_cc_id, back_face_attr.cc_id)){
          even_union_find.unify_sets(even_cc_id, back_face_attr.cc_id);
        }

        if (even_cc_id == nullptr && back_face_attr.cc_id != nullptr)
          even_cc_id = even_union_find.find(back_face_attr.cc_id);
      }
    }

    if (identified or beta3_identified){
      auto& face_attr = lcc.attribute<2>(face)->info();
      face_attr = DartInfo::FaceAttrValue();
      face_attr.plane[planeIteration] = true;
      face_attr.cc_id = odd_cc_id != nullptr ? odd_cc_id : odd_union_find.make_set(face);
    }

    if (identified) {
      assert(lcc.attribute<2>(back_face) == nullptr);
      assert(!lcc.is_free<3>(back_face));

      auto& back_face_attr = get_or_create_attr<2>(lcc, back_face)->info();
      back_face_attr.even = true;
      back_face_attr.plane[planeIteration] = true;
      back_face_attr.cc_id = even_cc_id != nullptr ? even_cc_id : even_union_find.make_set(lcc.beta<3>(back_face));
    }
  }

  void setup_next_level_plane(HexMeshingData& hdata){
    LCC& lcc = hdata.lcc;
    std::array<HexMeshingData::PlaneSet, 3> new_planes;

    for (int p = 0; p < 3; p++){
      HexMeshingData::PlaneSet new_plane_set;

      bool skip = false;
      int i = 0;
      for (auto& old_plane_set : hdata.first_face_of_planes[p]){
        if (!skip) { skip = true; continue;}

        std::queue<Dart_handle> to_explore;
        Union_find<Dart_handle> odd_union_find, even_union_find;

        size_type edge_mark = lcc.get_new_mark();
        size_type face_mark = lcc.get_new_mark();

        for (auto start : old_plane_set)
          to_explore.push(start);

        while (!to_explore.empty()){
          Dart_handle face = to_explore.front();
          to_explore.pop();

          setup_next_level_face(hdata, to_explore, odd_union_find, even_union_find, face, (Plane)p, edge_mark, face_mark);
        }

        HexMeshingData::PlaneCC odd_plane_cc, even_plane_cc;
        odd_plane_cc.reserve(odd_union_find.number_of_sets());
        even_plane_cc.reserve(odd_union_find.number_of_sets());

        for (auto uf_handle : get_partitions(odd_union_find)){
          Dart_handle face_cc = uf_handle->value;
          odd_plane_cc.push_back(face_cc);
        }

        for (auto uf_handle : get_partitions(even_union_find)){
          Dart_handle face_cc = uf_handle->value;
          even_plane_cc.push_back(face_cc);
        }

        new_plane_set.push_back(std::move(odd_plane_cc));
        new_plane_set.push_back(std::move(even_plane_cc));

        lcc.free_mark(edge_mark);
        lcc.free_mark(face_mark);
      }

      new_planes[p] = std::move(new_plane_set);
    }

    hdata.first_face_of_planes = std::move(new_planes);
  }

  void setup_next_level(HexMeshingData& hdata, MarkingFunction& cellIdentifier, Polyhedron& surface, Tree& aabb){
    setup_next_level_plane(hdata);
    clean_up_and_reevaluate_attributes(hdata, cellIdentifier, surface, aabb);
    hdata.level++;
  }

  void mark_identified_cells_from_3_attrs(HexMeshingData& hdata) {
    LCC& lcc = hdata.lcc;

    auto attributes = lcc.attributes<3>();

    for (auto it = attributes.begin(), end = attributes.end(); it != end; it++){
      if (it->info().type > static_cast<VolumeType>(0)){
        mark_all_0_cells<3>(lcc, it->dart(), hdata.identified_mark);
      }
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

    return Grid(pos, size, nb_subdiv_per_dim);
  }

  Grid generate_grid(LCC& lcc, Tree& aabb, int nb_subdiv_per_dim) {
    auto bbox = aabb.bbox();

    Point center = {bbox.xmin() + (bbox.x_span()/2),
                    bbox.ymin() + (bbox.y_span()/2),
                    bbox.zmin() + (bbox.z_span()/2)};

    double max_size = std::max(std::max(bbox.x_span(), bbox.y_span()), bbox.z_span());

    return generate_grid(lcc, center, max_size + 20, nb_subdiv_per_dim + 8);
  }

  void create_vertices_for_templates(HexMeshingData& hdata, RefinementData& rdata)
  {

    // 2 noeuds marqué l'un à coté de l'autre ne produit pas de sommet
    // 1 noeud marqué a coté d'un noeud non marqué produit un sommet

    // TODO A changer pour une itération sur les arretes ou pas selon comment l'algo va s'executer

    std::vector<Dart_handle> edges_to_subdivide;
    LCC& lcc = hdata.lcc;

    auto arrete_done = lcc.get_new_mark();

    int vertices_created = 0;
    for (auto dart : rdata.marked_nodes)
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

  bool mark_intersecting_volume_with_poly(LCC& lcc, Polyhedron& poly, Tree& tree, Dart_handle dart) {
    return is_intersect(lcc, dart, tree);
  }
}

namespace CGAL::HexRefinement {

  void render_two_refinement_result(const LCC& lcc, const char* title = "TwoRefinement Result"){
    LCCSceneOptions<LCC> gso;

    gso.draw_volume = [](const LCC& lcc, LCC::Dart_const_handle dart){return lcc.attribute<3>(dart) != nullptr;};
    gso.colored_volume = [](const LCC& lcc, LCC::Dart_const_handle dart){ return true; };
    gso.volume_color = [](const LCC& lcc, LCC::Dart_const_handle dart){
      auto& attr = lcc.attribute<3>(dart)->info();
      return attr.type == VolumeType::IDENTIFIED ? green() : blue();
    };

    CGAL::Graphics_scene buffer;
    add_to_graphics_scene(lcc, buffer, gso);
    CGAL::draw_graphics_scene(buffer);
  }

  void debug_render_2(const LCC& lcc, TwoRefinement::HexMeshingData& hdata, const char* title = "TwoRefinement Result"){
    LCCSceneOptions<LCC> gso;

    gso.draw_face = [](const LCC& lcc, LCC::Dart_const_handle dart){
      return true;
    };
    gso.face_color = [](const LCC& lcc, LCC::Dart_const_handle dart){
      auto a = lcc.attribute<2>(dart);
      auto b = lcc.attribute<3>(dart);
      return lcc.is_whole_cell_marked<2>(dart, debug3) ? red()
      : b != nullptr && b->info().type > VolumeType::REFINEMENT ? green()
      : blue();

      // return lcc.is_whole_cell_marked<3>(dart, debug3)? red() : a != nullptr && a->info().type == VolumeType::IDENTIFIED ? green() : blue();
    };
    gso.colored_face = [](const LCC& lcc, LCC::Dart_const_handle dart){ return true; };

    // gso.colored_vertex = [](const LCC& lcc, LCC::Dart_const_handle dart){
    //   return true;
    // };
    // gso.draw_vertex = [](const LCC& lcc, LCC::Dart_const_handle dart){
    //   return true;
    // };
    // gso.vertex_color = [&](const LCC& lcc, LCC::Dart_const_handle dart){
    //   return lcc.is_marked(dart, hdata.template_mark) ? red() : blue();
    // };
    CGAL::Graphics_scene buffer;
    add_to_graphics_scene(lcc, buffer, gso);
    CGAL::draw_graphics_scene(buffer);
  }

  void debug_render(const LCC& lcc, TwoRefinement::HexMeshingData& hdata, const char* title = "TwoRefinement Result"){
    LCCSceneOptions<LCC> gso;

    gso.draw_volume = [](const LCC& lcc, LCC::Dart_const_handle dart){
      return true;
    };
    gso.draw_face = [](const LCC& lcc, LCC::Dart_const_handle dart){ return true; };
    gso.face_color = [](const LCC& lcc, LCC::Dart_const_handle dart){
      auto a = lcc.attribute<2>(dart);
      auto b = lcc.attribute<3>(dart);
      return a != nullptr &&  a->info().even && a->info().plane[0]
      && b != nullptr && b->info().type > static_cast<VolumeType>(0) ? red()
      : b != nullptr && b->info().type > static_cast<VolumeType>(0) ? green() : blue();
      // return lcc.is_whole_cell_marked<2>(dart, debug3) ? red() :
      // b != nullptr && b->info().identified == 0 ? green() : blue();
      //: b != nullptr && b->info().identified == 0 ? green()
    };
    gso.colored_face = [](const LCC& lcc, LCC::Dart_const_handle dart){ return true; };

    CGAL::Graphics_scene buffer;
    add_to_graphics_scene(lcc, buffer, gso);
    CGAL::draw_graphics_scene(buffer);
  }

  LCC two_refinement(const std::string& file, int nb_subdiv_per_dim, TwoRefinement::MarkingFunction cellIdentifier, int nb_levels = 1) {
    using namespace TwoRefinement;

    HexMeshingData hdata;
    LCC& lcc = hdata.lcc;
    load_patterns(hdata.regular_templates, hdata.partial_templates);

    Polyhedron surface;
    Tree aabb;
    load_surface(file, surface, aabb);

    debug_node_mark = lcc.get_new_mark();
    debug_edge_mark = lcc.get_new_mark();
    debug3 = lcc.get_new_mark();

    hdata.grid = generate_grid(hdata.lcc, aabb, nb_subdiv_per_dim);
    hdata.identified_mark = lcc.get_new_mark();
    hdata.template_mark = lcc.get_new_mark();
    hdata.propagation_face_mark = lcc.get_new_mark();

    for (int i = 0; i < 2; i++){
      if (i == 0) setup_initial(hdata, cellIdentifier, surface, aabb);
      else setup_next_level(hdata, cellIdentifier, surface, aabb);

      expand_identified_cells(hdata, i, 2);

      for (int p = 0; p < 3; p++) {
        RefinementData rdata;

        mark_identified_cells_from_3_attrs(hdata);

        extract_darts_from_even_planes(hdata, rdata, (Plane)p);

        assert_dart_attr_are_unique<3>(lcc, rdata.volumes_to_refine, rdata.partial_templates_to_refine);

        create_vertices_for_templates(hdata, rdata);

        refine_marked_hexes(hdata, rdata);

        assert_all_faces_are_quadrilateral(lcc);
        assert_all_volumes_are_hexes(lcc);

        lcc.unmark_all(hdata.identified_mark);
        lcc.unmark_all(hdata.template_mark);
      }

      lcc.unmark_all(hdata.propagation_face_mark);

    }

    return lcc;
  }
}
