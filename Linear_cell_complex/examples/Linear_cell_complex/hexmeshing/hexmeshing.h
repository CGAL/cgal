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
#include <functional>

#include "../query_replace/cmap_query_replace.h"
#include "../query_replace/init_to_preserve_for_query_replace.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT       FT;
typedef Kernel::Point_3  Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Linear_cell_complex_for_combinatorial_map<3,3> LCC;
typedef typename LCC::Dart_handle Dart_handle;
typedef typename LCC::Vertex_attribute_handle Vertex_handle;
typedef typename LCC::size_type   size_type;

const size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

// to be removed
size_type debug_edge_mark;

/**
 * TODO: Make all functions templated for any LCC ?
 */

namespace CGAL::HexRefinement::TwoRefinement {
  enum Plane { XY, YZ, ZX, YX = XY, ZY = YZ, XZ = ZX };
  struct Grid { 
    Point center;
    float half_size;
    int nb_subdiv_per_dim;
  };

  std::vector<Dart_handle> marked_cells_in_hex(LCC& lcc, Dart_handle& dart, size_type mark){
    auto checked_mark = lcc.get_new_mark();
    auto darts_of_volume = lcc.darts_of_cell<3>(dart);
    std::vector<Dart_handle> marked;
    marked.reserve(4);

    for (auto it = darts_of_volume.begin(); it != darts_of_volume.end(); it++){
      Dart_handle dit = it;

      if (lcc.is_marked(dit, checked_mark) || !lcc.is_marked(dit, mark)) continue;

      lcc.mark_cell<0>(dit, checked_mark);
      marked.push_back(dit);
    };

    lcc.free_mark(checked_mark);
    return marked;
  }

  Dart_handle find_3_template_origin(LCC& lcc, std::vector<Dart_handle> &mcells_in_hex, size_type corner_mark) {
    Dart_handle d1 = mcells_in_hex[0];
    Dart_handle d2 = mcells_in_hex[1];
    Dart_handle d3 = mcells_in_hex[2];

    Dart_handle dart = d1;

    // If the two marks aren't on the same face, they must be one face away
    
    // TODO A améliorer, ex: itérerr sur chaque face 
    if (!lcc.belong_to_same_cell<2>(dart, d2)){
      dart = lcc.beta(dart, 2);
      assert(lcc.belong_to_same_cell<2>(dart, d2));
    }

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

  void refine_3_template(LCC &lcc, std::vector<Dart_handle> &mcells_in_hex, Pattern_substituer<LCC> &pss, size_type corner_mark)
  {
    Dart_handle origin_dart = find_3_template_origin(lcc, mcells_in_hex, corner_mark);  

    Dart_handle upper_d1 = origin_dart;
    Dart_handle upper_d2 = lcc.beta(origin_dart, 1, 1);

    assert(!lcc.is_marked(upper_d2, corner_mark));

    Dart_handle upper_edge = lcc.insert_cell_1_in_cell_2(upper_d1, upper_d2);

    // Query replace with the partial 3-template, making it into two volumes
    pss.query_replace_one_volume(lcc, origin_dart);

    // Face of the t<o neighboring volumes to the created volume
    Dart_handle face1 = lcc.beta(origin_dart, 2, 3);
    Dart_handle face2 = lcc.beta(origin_dart, 1, 2, 3);
    // lower_edge is created with the query_replace
    Dart_handle lower_edge = lcc.beta(upper_edge, 2, 1, 1);

    // Contract upper and lower edge into its barycenter 

    Dart_handle upper_mid_1 = lcc.insert_barycenter_in_cell<1>(upper_edge);
    Dart_handle upper_mid_2 = lcc.beta(upper_mid_1, 2, 1);
    Dart_handle lower_mid_1 = lcc.insert_barycenter_in_cell<1>(lower_edge);
    Dart_handle lower_mid_2 = lcc.beta(lower_mid_1, 2, 1);

    assert(!lcc.is_marked(upper_mid_2, corner_mark));
    assert(!lcc.is_marked(lower_mid_2, corner_mark));

    lcc.contract_cell<1>(upper_mid_1);
    lcc.contract_cell<1>(upper_mid_2);
    lcc.contract_cell<1>(lower_mid_1);
    lcc.contract_cell<1>(lower_mid_2);

    // Contract the two remaining 2-darts faces 

    Dart_handle face_to_remove_1 = lcc.beta(upper_edge, 2, 1);
    Dart_handle face_to_remove_2 = lcc.beta(upper_edge, 2, 1, 1, 2, 1, 2);

    assert(lcc.darts_of_orbit<1>(face_to_remove_1).size() == 2);
    assert(lcc.darts_of_orbit<1>(face_to_remove_2).size() == 2);

    lcc.contract_cell<2>(face_to_remove_1);
    lcc.contract_cell<2>(face_to_remove_2);

    // Remove the created volume and sew the two neighboring volumes
    lcc.remove_cell<3>(upper_d1);
    lcc.sew<3>(face1, face2);
  }

  void refine_marked_hexes(LCC &lcc, Pattern_substituer<LCC> &ps_regular, Pattern_substituer<LCC> &ps_partial_3t, size_type template_mark, size_type corner_mark)
  {
    // Replace all faces that has a matching fpattern
    std::vector<Dart_handle> marked_cells;


    for (auto dit = lcc.one_dart_per_cell<2>().begin(),
              dend = lcc.one_dart_per_cell<2>().end();
        dit != dend;
        dit++)
    {
      Dart_handle dart = dit;
      marked_cells.push_back(dart);
    }


    int nbsub = 0;
    for (auto& dart : marked_cells)
    {
      if (ps_regular.query_replace_one_face(lcc, dart, template_mark) != SIZE_T_MAX) nbsub++;
    }

    std::cout << nbsub << " Faces substitution was made" << std::endl;

    // Refine the resulting volumes
    marked_cells.clear();

    using Volume_3_Set = std::vector<std::vector<Dart_handle>>;
    std::vector<Dart_handle>& volumes_regular_template = marked_cells;
    Volume_3_Set volumes_3_template;

    for (auto dit = lcc.one_dart_per_cell<3>().begin(),
              dend = lcc.one_dart_per_cell<3>().end();
        dit != dend;
        dit++)
    {
      Dart_handle dart = dit;
      std::vector mcell_in_hex = marked_cells_in_hex(lcc, dart, template_mark);

      if (mcell_in_hex.size() == 3){
        volumes_3_template.push_back(mcell_in_hex);
        continue;
      }

      marked_cells.push_back(dart);
    }

    nbsub = 0;
    for (auto& dart : volumes_regular_template)
    {
      if (ps_regular.query_replace_one_volume(lcc, dart) != SIZE_T_MAX) nbsub++;
    }

    // Refine remaining 3 patterns

    // Find the four middle vertices
    for (auto mcells_in_hex : volumes_3_template){
      refine_3_template(lcc, mcells_in_hex, ps_partial_3t, corner_mark);
      nbsub++;
    }

    std::cout << nbsub << " volumic substitution was made" << std::endl;
  }

  void __mark_line(LCC& lcc, Dart_handle dit, std::vector<Dart_handle> &marked_cells, size_type template_mark, int nb_subdiv){
    // Mark the first cell
    lcc.mark_cell<0>(dit, template_mark);

    // Mark the line 
    for (int i = 0; i < nb_subdiv; i++) {
      // Dart on the next vertex
      dit = lcc.beta(dit, 0);
      
      // Change volume (if exists)
      if (i < nb_subdiv - 1) 
        dit = lcc.beta(dit, 0, 2, 3, 2); 
      
      lcc.mark_cell<0>(dit, template_mark);
      marked_cells.push_back(dit);
    }
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

  // renommer  : between cell pairs
  std::vector<Dart_handle> mark_nodes_in_cell_pair(LCC& lcc, Grid grid, Plane plane, size_type identified_mark, size_type template_mark){
    std::vector<Dart_handle> marked_cells;
    
    // Method:  iterate line by line the studied plane (2d), 
    // To get the next vertex on the line, we use the third dimension, because it is garanted to have an other plane of volume behind
    // because we are operating on every even plane of the grid. 

    // the first dart is the lowest in coordinates (-x, -y, -z), pointing towards -z; on the ZX plane  
    Dart_handle start_plane = __first_vertex_of_even_plane(lcc, plane);
    Dart_handle start_line = start_plane; // plan xy

    int nb_subdiv = grid.nb_subdiv_per_dim;

    for (int z = 0; z < nb_subdiv / 2; z++) {
      for (int x = 0; x < nb_subdiv; x++){
        __mark_line(lcc, start_line, marked_cells, template_mark, nb_subdiv);
        
        // Walk to the next line 
        start_line = lcc.beta(start_line, 2, 0, 0, 2); 

        if (x < nb_subdiv - 1) start_line = lcc.beta(start_line, 3);
      }
      
      // The last starting dart is of opposite direction, walk backward into the odd plane
      // for the dart to be incident from the proper vertex
      start_line = lcc.beta(start_line, 1, 2, 3, 2, 1);
      __mark_line(lcc, start_line, marked_cells, template_mark, nb_subdiv);

      // Walk to the next even plane
      start_plane = __next_even_plane(lcc, start_plane, plane); 
      start_line = start_plane;
    }

    return marked_cells;

  }

  Grid generate_grid(LCC& lcc, Point center, float half_size, int nb_subdiv_per_dim) {
    float s = (2*half_size) / nb_subdiv_per_dim;

    CGAL_precondition_msg((nb_subdiv_per_dim % 2) == 0, "nb_subdiv must be even");

    // hexahedron translated by pos
    Vector pos;

    for (int x = 0; x < nb_subdiv_per_dim; x++) {
      for (int y = 0; y < nb_subdiv_per_dim; y++) {
        for (int z = 0; z < nb_subdiv_per_dim; z++) {
          pos = {x * s, y * s, z * s};
          lcc.make_hexahedron(Point(0,0,0) + pos, Point(s,0,0) + pos,
                              Point(s,s,0) + pos, Point(0,s,0) + pos,
                              Point(0,s,s) + pos, Point(0,0,s) + pos,
                              Point(s,0,s) + pos, Point(s,s,s) + pos);
        }
      }
    }

    lcc.sew3_same_facets();

    return Grid{.center = center, .half_size = half_size, .nb_subdiv_per_dim = nb_subdiv_per_dim };
  }

  void create_vertices_for_templates(LCC &lcc, std::vector<Dart_handle>& marked_cells, size_type template_mark)
  {

    // 2 noeuds marqué l'un à coté de l'autre ne produit pas de sommet
    // 1 noeud marqué a coté d'un noeud non marqué produit un sommet

    // TODO A changer pour une itération sur les arretes ou pas selon comment l'algo va s'executer

    std::vector<Dart_handle> edges_to_subdivide;

    auto arrete_done = lcc.get_new_mark();

    for (auto dart : marked_cells)
    {
      for (auto nit = lcc.one_dart_per_incident_cell<1, 0>(dart).begin(),
                nend = lcc.one_dart_per_incident_cell<1, 0>(dart).end();
          nit != nend;
          nit++)
      {
        if (lcc.is_marked(nit, arrete_done))
          continue;

        // If the node is next to an other marked node, we don't have to create vertices
        if (lcc.is_marked(lcc.beta<1>(nit), template_mark)){
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

  void load_patterns(Pattern_substituer<LCC> &regular_templates, Pattern_substituer<LCC>& partial_3_template) {
    // BUG Manque des opérateurs de déplacements pour Pattern, le reserve est un fix temporaire
    // Pour pouvoir charger les patterns correctement sans réallocation
    regular_templates.m_fpatterns.reserve(10);
    regular_templates.load_additional_fpattern(CGAL::data_file_path("hexmeshing/fpattern/pattern1-face.moka"), mark_1template_face);
    regular_templates.load_additional_fpattern(CGAL::data_file_path("hexmeshing/fpattern/pattern2-face.moka"), mark_2template_face);
    regular_templates.load_vpatterns(CGAL::data_file_path("hexmeshing/vpattern/complete"));
  
    // TODO: Chargé séparément pour le moment, voir si je laisse comme ça 
    partial_3_template.load_vpatterns(CGAL::data_file_path("hexmeshing/vpattern/partial"));
  }
} 

namespace CGAL::HexRefinement {

  // Identifies which 3-cell should be refined 
  using MarkingFunction = std::function<void(LCC&, size_type)>;

  void tworefinement(LCC& lcc, MarkingFunction cellIdentifier) {
    using namespace TwoRefinement;

    LCC lcc_grid;
    
    Pattern_substituer<LCC> ps_regular, ps_partial_3t;
    load_patterns(ps_regular, ps_partial_3t);

  
    auto identified_mark = lcc.get_new_mark();
    auto template_mark = lcc.get_new_mark();
    auto corner_mark = lcc.get_new_mark();
    lcc.negate_mark(corner_mark);


    Grid grid = generate_grid(lcc_grid, Point(0,0,0), 5, 10);

    cellIdentifier(lcc_grid, identified_mark);

    std::vector<Dart_handle> marked_cells = mark_nodes_in_cell_pair(lcc, grid, Plane::XY, identified_mark, template_mark);
    create_vertices_for_templates(lcc, marked_cells, template_mark);
    refine_marked_hexes(lcc, ps_regular, ps_partial_3t, template_mark, corner_mark);

    lcc.free_mark(corner_mark);
  } 

  // temporaire
  void tworefinement_one_step_from_marked(LCC& lcc, std::vector<Dart_handle> marked_cells, size_type template_mark) {
    using namespace TwoRefinement;

    Pattern_substituer<LCC> ps_regular, ps_partial_3t;
    load_patterns(ps_regular, ps_partial_3t);
    
    auto corner_mark = lcc.get_new_mark();
    lcc.negate_mark(corner_mark);

    create_vertices_for_templates(lcc, marked_cells, template_mark);
    refine_marked_hexes(lcc, ps_regular, ps_partial_3t, template_mark, corner_mark);

    lcc.free_mark(corner_mark);
  }
}
