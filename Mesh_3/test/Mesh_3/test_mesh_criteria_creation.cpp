#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>

// Domain
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;

// Triangulation
typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
// Criteria
typedef CGAL::Mesh_criteria_3<Tr> Mc;
typedef CGAL::Mesh_constant_domain_field_3<Tr::Geom_traits,
                                           Mesh_domain::Index> Sizing_field;

typedef Sizing_field  Esf;
typedef Sizing_field  Fsf;
typedef Sizing_field  Csf;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

int main()
{
  Tr::Point p1(0,0,0);
  Tr::Point p2(1,0,0);
  Tr::Point p3(0,1,0);
  Tr::Point p4(0,0,1);

  Tr tr;
  tr.insert(p1);
  tr.insert(p2);
  tr.insert(p3);
  tr.insert(p4);

  Tr::Cell_handle ch = tr.finite_cells_begin();
  int k = 0;
  Tr::Facet f = std::make_pair(ch,k);
  Tr::Cell::Subdomain_index sub_index = 1;
  Tr::Cell::Surface_patch_index surf_index = 2;
  Tr::Cell::Surface_patch_index surf_index_bis = 21;
  Tr::Vertex::Index index (surf_index);
  Tr::Vertex::Index index_bis(surf_index_bis);

  // Init cell
  tr.dual(ch);
  ch->set_subdomain_index(sub_index);

  // Init facet
  Tr::Bare_point facet_circum =
    tr.geom_traits().construct_weighted_circumcenter_3_object()(
      tr.point(ch, k+1), tr.point(ch, k+2), tr.point(ch, k+3));

  ch->set_surface_patch_index(k,surf_index);
  ch->set_facet_surface_center(k,facet_circum);
  ch->set_facet_surface_center_index(k,index);

  // Init vertices
  ch->vertex(0)->set_dimension(2);
  ch->vertex(1)->set_dimension(2);
  ch->vertex(2)->set_dimension(2);
  ch->vertex(3)->set_dimension(2);

  ch->vertex(0)->set_index(index);
  ch->vertex(1)->set_index(index);
  ch->vertex(2)->set_index(index_bis);
  ch->vertex(3)->set_index(index_bis);

  // -----------------------------------
  // Test edge criteria
  // -----------------------------------
  Tr::Bare_point bp1 = tr.geom_traits().construct_point_3_object()(p1);

  Mc ec1(edge_size = 1);
  assert( ec1.edge_criteria_object().sizing_field(bp1,1,index) == 1 );

  Mc ec2(edge_size = Esf(2));
  assert( ec2.edge_criteria_object().sizing_field(bp1,1,index) == 2 );

  Mc ec3(edge_size = 3.);
  assert( ec3.edge_criteria_object().sizing_field(bp1,1,index) == 3 );

  Mc ec8(edge_distance = 8.);
  assert( ec8.edge_criteria_object().distance_field(bp1, 1, index) == 8. );

  Mc ec9(edge_distance = Esf(9.));
  assert( ec9.edge_criteria_object().distance_field(bp1, 1, index) == 9.);


  // -----------------------------------
  // Test facet criteria
  // -----------------------------------
  typedef Tr::Geom_traits::FT FT;
  Tr::Geom_traits::Compute_squared_radius_3 squared_radius = tr.geom_traits().compute_squared_radius_3_object();
  Tr::Geom_traits::Construct_point_3 cp = tr.geom_traits().construct_point_3_object();

  FT radius_facet = CGAL::sqrt(squared_radius(cp(tr.point(ch, k+1)),
                                              cp(tr.point(ch, k+2)),
                                              cp(tr.point(ch, k+3))));

  FT facet_size_ok = radius_facet*FT(10);

  Mc fc1(facet_size = facet_size_ok);
  assert( ! fc1.facet_criteria_object()(tr, f) );

  Mc fc3(facet_size = Fsf(facet_size_ok));
  assert( ! fc3.facet_criteria_object()(tr, f) );

  Mc fc8(facet_distance = 8.);
  Mc fc9(facet_angle = 9.);
  Mc fc10(facet_angle = 10.1,
          facet_distance = 10.2,
          facet_size = 10.3,
          facet_min_size = 0.2);
  Mc fc1Ob(facet_angle = 10.11,
          facet_distance = 10.12,
          facet_min_size = 0.3,
          facet_size = Fsf(10.14));

  // Test construction from int
  Mc fc11(facet_size = 11);
  Mc fc11b(facet_size = 11, facet_min_size = 1);
  Mc fc12(facet_size = Fsf(12));
  Mc fc12b(facet_size = Fsf(12), facet_min_size = 1);

  // Test topological criterion creation
  Mc fc14(facet_topology = CGAL::FACET_VERTICES_ON_SURFACE);
  assert( ! fc14.facet_criteria_object()(tr, f) );

  Mc fc15(facet_topology = CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH);
  assert( fc15.facet_criteria_object()(tr, f) );

  // -----------------------------------
  // Test cell criteria
  // -----------------------------------
  FT radius_cell = CGAL::sqrt(squared_radius(cp(tr.point(ch, 0)),
                                             cp(tr.point(ch, 1)),
                                             cp(tr.point(ch, 2)),
                                             cp(tr.point(ch, 3))));

  FT cell_size_ok = radius_cell*FT(10);

  Mc cc1(cell_size = cell_size_ok);
  assert( ! cc1.cell_criteria_object()(tr, ch) );

  Mc cc3(cell_size = Fsf(cell_size_ok));
  assert( ! cc3.cell_criteria_object()(tr, ch) );

  Mc cc7(cell_size = Csf(cell_size_ok));
  assert( ! cc7.cell_criteria_object()(tr, ch) );

  Mc cc8(cell_radius_edge_ratio = 8.);
  Mc cc9(cell_radius_edge_ratio = 9.1,
         cell_size = Csf(9.2) );
  Mc cc10(cell_radius_edge_ratio = 10.1,
          cell_size = 10.2,
          cell_min_size = 0.1);
  Mc cc10b(cell_radius_edge_ratio = 10.1,
          cell_min_size = 0.1,
          cell_size = Csf(10.3));

  // Test construction from int
  Mc cc11(cell_size = 11);
  Mc cc11b(cell_size = 11, cell_min_size = 1);
}
