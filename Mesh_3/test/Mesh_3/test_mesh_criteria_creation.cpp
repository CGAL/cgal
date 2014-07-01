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
  Tr::Point facet_circum =
    tr.geom_traits().construct_weighted_circumcenter_3_object()(
      ch->vertex(k+1)->point(), ch->vertex(k+2)->point(), ch->vertex(k+3)->point());
  
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
  Mc ec1(edge_size = 1);
  assert( ec1.edge_criteria_object().sizing_field(p1,1,index) == 1 );
  
  Mc ec2(edge_sizing_field = Esf(2));
  assert( ec2.edge_criteria_object().sizing_field(p1,1,index) == 2 );

  Mc ec3(edge_sizing_field = 3.);
  assert( ec3.edge_criteria_object().sizing_field(p1,1,index) == 3 );
  
  Mc ec4(edge_size = 4.1,
         edge_sizing_field = Esf(4.2));
  assert( ec4.edge_criteria_object().sizing_field(p1,1,index) == 4.1 );
  
  Mc ec5(sizing_field = 5.);
  assert( ec5.edge_criteria_object().sizing_field(p1,1,index) == 5 );
  
  Mc ec6(sizing_field = 6.1,
         edge_sizing_field = 6.2);
  assert( ec6.edge_criteria_object().sizing_field(p1,1,index) == 6.2 );
  
  Mc ec7(sizing_field = 7.1,
         edge_size = 7.2);
  assert( ec7.edge_criteria_object().sizing_field(p1,1,index) == 7.2 );
  
  
  // -----------------------------------
  // Test facet criteria
  // -----------------------------------
  typedef Tr::Geom_traits::FT FT;
  FT radius_facet = CGAL::sqrt(CGAL::squared_radius(ch->vertex(k+1)->point(),
                                                    ch->vertex(k+2)->point(),
                                                    ch->vertex(k+3)->point()));
  
  FT facet_size_ok = radius_facet*FT(10);
  FT facet_size_nok = radius_facet/FT(10);

  Mc fc1(facet_size = facet_size_ok);
  assert( ! fc1.facet_criteria_object()(f) );
  
  Mc fc2(facet_sizing_field = facet_size_ok);
  assert( ! fc2.facet_criteria_object()(f) );

  Mc fc3(facet_sizing_field = Fsf(facet_size_ok));
  assert( ! fc3.facet_criteria_object()(f) );

  Mc fc4(facet_sizing_field = facet_size_nok,
         facet_size = facet_size_ok);
  assert( ! fc4.facet_criteria_object()(f) );
  
  Mc fc5(sizing_field = facet_size_ok);
  assert( ! fc5.facet_criteria_object()(f) );
  
  Mc fc6(facet_size = facet_size_ok,
         facet_sizing_field = facet_size_nok,
         sizing_field = facet_size_nok);
  assert( ! fc6.facet_criteria_object()(f) );
  
  Mc fc7(facet_sizing_field = Fsf(facet_size_ok),
         sizing_field = facet_size_nok);
  assert( ! fc7.facet_criteria_object()(f) );
  
  Mc fc8(facet_distance = 8.);
  Mc fc9(facet_angle = 9.);
  Mc fc10(facet_angle = 10.1,
          facet_distance = 10.2,
          facet_size = 10.3,
          facet_sizing_field = Fsf(10.4),
          sizing_field = 10.5);
  
  // Test construction from int
  Mc fc11(facet_size = 11);
  Mc fc12(facet_sizing_field = 12);
  Mc fc13(sizing_field = 13);
  
  // Test topological criterion creation
  Mc fc14(facet_topology = CGAL::FACET_VERTICES_ON_SURFACE);
  assert( ! fc14.facet_criteria_object()(f) );
  
  Mc fc15(facet_topology = CGAL::FACET_VERTICES_ON_SAME_SURFACE_PATCH);
  assert( fc15.facet_criteria_object()(f) );
  
  // -----------------------------------
  // Test cell criteria
  // -----------------------------------
  FT radius_cell = CGAL::sqrt(CGAL::squared_radius(ch->vertex(0)->point(),
                                                   ch->vertex(1)->point(),
                                                   ch->vertex(2)->point(),
                                                   ch->vertex(3)->point()));
  
  FT cell_size_ok = radius_cell*FT(10);
  FT cell_size_nok = radius_cell/FT(10);
  
  Mc cc1(cell_size = cell_size_ok);
  assert( ! cc1.cell_criteria_object()(ch) );
  
  Mc cc2(cell_sizing_field = cell_size_ok);
  assert( ! cc2.cell_criteria_object()(ch) );
  
  Mc cc3(cell_sizing_field = Fsf(cell_size_ok));
  assert( ! cc3.cell_criteria_object()(ch) );
  
  Mc cc4(cell_sizing_field = cell_size_nok,
         cell_size = cell_size_ok);
  assert( ! cc4.cell_criteria_object()(ch) );
  
  Mc cc5(sizing_field = cell_size_ok);
  assert( ! cc5.cell_criteria_object()(ch) );
  
  Mc cc6(cell_size = cell_size_ok,
         cell_sizing_field = cell_size_nok,
         sizing_field = cell_size_nok);
  assert( ! cc6.cell_criteria_object()(ch) );
  
  Mc cc7(cell_sizing_field = Csf(cell_size_ok),
         sizing_field = cell_size_nok);
  assert( ! cc7.cell_criteria_object()(ch) );
  
  Mc cc8(cell_radius_edge_ratio = 8.);
  Mc cc9(cell_radius_edge_ratio = 9.1,
         sizing_field = Csf(9.2) );
  Mc cc10(cell_radius_edge_ratio = 10.1,
          cell_size = 10.2,
          cell_sizing_field = Csf(10.3),
          sizing_field = 10.4);
  
  // Test construction from int
  Mc cc11(cell_size = 11);
  Mc cc12(cell_sizing_field = 12);
  Mc cc13(sizing_field = 13);
}
