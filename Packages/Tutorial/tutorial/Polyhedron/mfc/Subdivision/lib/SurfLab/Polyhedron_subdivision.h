// ======================================================================
//
// Copyright (c) 2002 SurfLab of CISE of University of Florida
//
// File          : libs/src/cgalExt/Polyhedron_subdivision.h
// Description   : Provides auxiliary functions to subdivide a polyhedron
// Creation_date : 29 Jan 2002
// Author(s)     : Le-Jeng Shiue <sle-jeng@cise.ufl.edu>
//
// ======================================================================

// $Id: Polyhedron_subdivision.h,v 1.2 2002/01/25 19:01:22 sle-jeng Exp 

#ifndef _POLYHEDRON_SUBDIVISION_H_01292002
#define _POLYHEDRON_SUBDIVISION_H_01292002

#include <vector>

#include <CGAL/circulator.h>

#include "Polyhedron_decorator.h"
#include "Polyhedron_subdivision_rules.h"
#include "Polyhedron_memory_builder.h"

// ======================================================================
///
template <class _Poly>
class Polyhedron_subdivision : public Polyhedron_decorator<_Poly> {
  typedef _Poly                                        Polyhedron;

  typedef typename Polyhedron::Traits                  Traits;
  typedef typename Traits::Kernel                      Kernel;

  typedef typename Polyhedron::Halfedge_data_structure HDS;
  typedef typename Polyhedron::Vertex                  Vertex;
  typedef typename Polyhedron::Halfedge                Halfedge;
  typedef typename Polyhedron::Facet                   Facet;
  
  typedef typename Polyhedron::Vertex_handle           Vertex_handle;
  typedef typename Polyhedron::Halfedge_handle         Halfedge_handle;
  typedef typename Polyhedron::Facet_handle            Facet_handle;

  typedef typename Polyhedron::Vertex_iterator         Vertex_iterator;
  typedef typename Polyhedron::Halfedge_iterator       Halfedge_iterator;
  typedef typename Polyhedron::Edge_iterator           Edge_iterator;
  typedef typename Polyhedron::Facet_iterator          Facet_iterator;

  typedef typename Polyhedron::Halfedge_around_facet_circulator  
                                            Halfedge_around_facet_circulator;
  typedef typename Polyhedron::Halfedge_around_vertex_circulator 
                                            Halfedge_around_vertex_circulator;

  typedef typename Polyhedron::Point_3                 Point;
  typedef typename Kernel::FT                          FT;

/**@name Class Methods */
//@{
public:
  /** Subdivide the input polyhedron p with Catmull-Clark subdivision 
      rules.
      
      Precondition: 
      Postcondition:
      The input  polyhedron p is subdivided to the limit surface of 
      Catmull-Clark subdivision.
  */
  static void CatmullClark_subdivision(Polyhedron& p, int step = 1) {
    quad_quadralize_polyhedron(p, CatmullClark_rule<Polyhedron>(), step);
  }

  /** Subdivide the input polyhedron p with Loop subdivision 
      rules.
  */
  static void Loop_subdivision(Polyhedron& p, int step = 1) {
    tri_quadralize_polyhedron(p, Loop_rule<Polyhedron>() , step);
  }
  /** Subdivide the input polyhedron p with Loop subdivision 
      rules.
  */
  static void QT43_subdivision(Polyhedron& p, int step = 1) {
    qt_quadralize_polyhedron(p, QT43_rule<Polyhedron>(), step);
  }

  ///
  static void PTQ(Polyhedron& p, int step = 1) {
    tri_quadralize_polyhedron(p, average_rule<Polyhedron>(), step);
  }
  ///
  static void PQQ(Polyhedron& p, int step = 1) {
    quad_quadralize_polyhedron(p, average_rule<Polyhedron>(), step);
  }
  ///
  static void PQTQ(Polyhedron& p, int step = 1) {
    qt_quadralize_polyhedron(p, average_rule<Polyhedron>(), step);
  }

  /** Make all faces in the input polyhedron p valance 4 with user's
      specialized rule. MT indicates the type of the refined mesh: 
      3(tris) or 4 (quads)
      
      Precondition: 
      Postcondition:
      The input  polyhedron p has only regular face. All new 
      edge-vertices are the average of their edge-neighbor-vertices.
      All new face-vertices are the average of their face-neighbor-
      vertices.
  */
//   template <class RULE, int MT>
//   static void quadralize_polyhedron(Polyhedron& p, int step = 1) {
//     // Check the template classes are matching to each other
//     typedef typename RULE::Polyhedron RPolyhedron;
//     assert_equal_types(RPolyhedron(), Polyhedron());

//     if (MT == 4) for (int i = 0; i < step; i++) quad_quadralize_1step<RULE>(p);
//     if (MT == 3) for (int i = 0; i < step; i++) tri_quadralize_1step<RULE>(p);
//     if (MT == 43) for (int i = 0; i < step; i++) qt_quadralize_1step<RULE>(p);
//   }
  ///
  template <template <typename> class RULE>
  static void quad_quadralize_polyhedron(Polyhedron& p, RULE<Polyhedron> rule, int step = 1) {
    for (int i = 0; i < step; i++) quad_quadralize_1step(p, rule);
  }
  ///
  template <template <typename> class RULE>
  static void tri_quadralize_polyhedron(Polyhedron& p, RULE<Polyhedron> rule, int step = 1) {
    for (int i = 0; i < step; i++) tri_quadralize_1step(p, rule);
  }

  ///
  template <template <typename> class RULE>
  static void qt_quadralize_polyhedron(Polyhedron& p, RULE<Polyhedron> rule, int step = 1) {
    for (int i = 0; i < step; i++) qt_quadralize_1step(p, rule);
  }


  /** Subdivide the input polyhedron p with Doo-Sabin subdivision 
      rules.
      
      Precondition: 
      Postcondition:
      The input  polyhedron p is subdivided to the limit surface of 
      Doo-Sabin subdivision.
  */
  static void DooSabin_subdivision(Polyhedron& p, int step = 1) {
    dualize_polyhedron(p, DooSabin_rule<Polyhedron>(), step);
  }

  /** Dualize the input polyhedron p. Make all vertex in the input 
      polyhedron p valance 4 with user's specialized rule. 
      
      Precondition: 
      Postcondition:
  */
  template <template <typename> class RULE>
  static void dualize_polyhedron(Polyhedron& p, RULE<Polyhedron> rule, int step = 1) {
    for (int i = 0; i < step; ++i) dualize_1step(p, rule);
  }


  /** Subdivide the input polyhedron p with Sqrt(3) subdivision 
      rules.
      
      Precondition: 
      Postcondition:
      The input  polyhedron p is subdivided to the limit surface of 
      Sqrt(3) subdivision.
  */
  static void Sqrt3_subdivision(Polyhedron& p, int step = 1) {
    sqrt3refine_polyhedron(p, Sqrt3_rule<Polyhedron>(), step);
  }

  /** Make all faces in the input polyhedron p valance 4 with user's
      specialized rule. MT indicates the type of the refined mesh: 
      3(tris) or 4 (quads)
      
      Precondition: 
      Postcondition:
      The input  polyhedron p has only regular face. All new 
      edge-vertices are the average of their edge-neighbor-vertices.
      All new face-vertices are the average of their face-neighbor-
      vertices.
  */
  template <template <typename> class RULE>
  static void sqrt3refine_polyhedron(Polyhedron& p, RULE<Polyhedron> rule, int step = 1) {
    for (int i = 0; i < step; i++) sqrt3refine_1step(p, rule);
  }



protected:
  /// It actually can take care non-quad facets
  template <template <typename> class RULE>
  static void quad_quadralize_1step(Polyhedron& p, RULE<Polyhedron> rule);
  ///
  template <template <typename> class RULE>
  static void tri_quadralize_1step(Polyhedron& p, RULE<Polyhedron> rule);
  ///
  template <template <typename> class RULE>
  static void qt_quadralize_1step(Polyhedron& p, RULE<Polyhedron> rule);

  ///
  template <template <typename> class RULE>
  static void dualize_1step(Polyhedron& p, RULE<Polyhedron> rule);

  ///
  template <template <typename> class RULE>
  static void sqrt3refine_1step(Polyhedron& p, RULE<Polyhedron> rule);

//@}
};


// ======================================================================
///
template <class _P> template <template <typename> class RULE>
void Polyhedron_subdivision<_P>::quad_quadralize_1step(_P& p, RULE<_P> rule) {
  p.normalize_border();

  // Build a new vertices buffer has the following structure
  //
  // 0 1 ... e_begin ... f_begin ... (end_of_buffer)
  // 0 ... e_begin-1       : store the positions of the vertex-vertices
  // e_begin ... f_begin-1 : store the positions of the edge-vertices
  // f_begin ... (end)     : store the positions of the face-vertices
  // The index of the vertices buffer should 1-1 map to the distance
  // of the corresponding iterator to the begin of the iterator.
  int num_vertex = p.size_of_vertices();
  int num_edge = p.size_of_halfedges()/2;
  int num_facet = p.size_of_facets();

  // If Polyhedron is using vector, we need to reserve the memory to prevent 
  // the CGAL_assertion.
  // This function for polyhedron using list is VOID.
  p.reserve(num_vertex+num_edge+num_facet, 4*2*num_edge, 4*num_edge/2);

  Point* vertex_point_buffer = new Point[num_vertex + num_edge + num_facet];
  Point* edge_point_buffer = vertex_point_buffer + num_vertex;
  Point* face_point_buffer = edge_point_buffer + num_edge;

  std::vector<bool> v_onborder(num_vertex);

  Facet_iterator fitr = p.facets_begin();
  for (int i = 0; i < num_facet; i++, ++fitr)
    rule.face_point_rule(fitr, face_point_buffer[i]);

  int sb = p.size_of_border_edges();

  Edge_iterator eitr = p.edges_begin();
  for (int i = 0; i < num_edge-sb; i++, ++eitr)
    rule.edge_point_rule(eitr, edge_point_buffer[i]);
  for (int i = num_edge-sb; i < num_edge; i++, ++eitr) {
    int v = std::distance(p.vertices_begin(), eitr->vertex());
    v_onborder[v] = true;
    rule.border_point_rule(eitr, edge_point_buffer[i], vertex_point_buffer[v]);
  }

  Vertex_iterator vitr = p.vertices_begin();
  for (int i = 0; i < num_vertex; i++, ++vitr)
    if (!v_onborder[i]) rule.vertex_point_rule(vitr, vertex_point_buffer[i]);

  // Build the connectifty using insert_vertex() and insert_edge()
  // 1. insert_vertex() to all edges and set them to new positions
  // 2. insert_edge() between 2 randomly selected neighboring new inserted 
  //    vertices
  // 3. insert_vertex() to the new inserted edge and set them to new positions
  // 4. insert_edge() between all other new inserted vetices of step 1 and
  //    the new inserted vertex of step 3
  // Step 1.
  eitr = p.edges_begin();
  for (int i = 0; i < num_edge; i++, ++eitr) {
    Vertex_handle vh = insert_vertex(p, eitr);
    vh->point() = edge_point_buffer[i];
  }
  fitr = p.facets_begin();

  // TODO: the topoloy modification can be done by a template function
  //       and that gives the user a chance to create new topoly rules.
  for (int i = 0; i < num_facet; i++, ++fitr) {
    // Step 2.
    Halfedge_around_facet_circulator hcir_begin = fitr->facet_begin();
    Halfedge_around_facet_circulator hcir = hcir_begin;
    
    Halfedge_handle e1 = ++hcir; // e1 points to the newly inserted vertex
    ++hcir; // Skips one original vertex 
    Halfedge_handle e2 = ++hcir; // points to the next newly inserted vertex
    ++hcir; // Must move the cir before inserts the new edge !!
    Halfedge_handle newe = insert_edge(p, e1, e2);
    
    // Step 3.
    Halfedge_handle newv = insert_vertex_return_edge(p, newe);
    newv = newv->opposite()->prev(); // change newv to the larger face and 
                                     // still points to the newly inserted 
                                     // vertex
    // Update the geometry data of the newly inserted face-vertices
    newv->vertex()->point() = face_point_buffer[i];
   
    // Step 4.
    while (hcir != hcir_begin) {
      e1 = ++hcir;
      ++hcir; // Must move the cir before inserts the new edge !!
      insert_edge(p, e1, newv); 
    }
  }
  
  // Update the geometry data of the newly inserted vertices by the 
  // vertices buffer
  vitr = p.vertices_begin();
  for (int i = 0; i < num_vertex; i++, ++vitr) 
    vitr->point() = vertex_point_buffer[i];

  delete []vertex_point_buffer;
}

// ======================================================================
///
template <class _P> template <template <typename> class RULE>
void Polyhedron_subdivision<_P>::tri_quadralize_1step(_P& p, RULE<_P> rule) {
  p.normalize_border();

  // Build a new vertices buffer has the following structure
  //
  // 0 1 ... e_begin ... f_begin ... (end_of_buffer)
  // 0 ... e_begin-1       : store the positions of the vertex-vertices
  // e_begin ... (end)     : store the positions of the edge-vertices
  // The index of the vertices buffer should 1-1 map to the distance
  // of the corresponding iterator to the begin of the iterator.
  int num_vertex = p.size_of_vertices();
  int num_edge = p.size_of_halfedges()/2;
  int num_facet = p.size_of_facets();

  // If Polyhedron is using vector, we need to reserve the memory to prevent 
  // the CGAL_assertion.
  // This function for polyhedron using list is VOID.
  p.reserve(num_vertex+num_edge, 2*2*num_edge, 4*num_edge/2);

  Point* vertex_point_buffer = new Point[num_vertex + num_edge];
  Point* edge_point_buffer = vertex_point_buffer + num_vertex;

  std::vector<bool> v_onborder(num_vertex);
  int sb = p.size_of_border_edges();

  Edge_iterator eitr = p.edges_begin();
  for (int i = 0; i < num_edge-sb; i++, ++eitr)
    rule.edge_point_rule(eitr, edge_point_buffer[i]);
  for (int i = num_edge-sb; i < num_edge; i++, ++eitr) {
    int v = std::distance(p.vertices_begin(), eitr->vertex());
    v_onborder[v] = true;
    rule.border_point_rule(eitr, edge_point_buffer[i], vertex_point_buffer[v]);
  }

  Vertex_iterator vitr = p.vertices_begin();
  for (int i = 0; i < num_vertex; i++, ++vitr)
    if (!v_onborder[i]) rule.vertex_point_rule(vitr, vertex_point_buffer[i]);

  // Build the connectifty using insert_vertex() and insert_edge()
  // 1. insert_vertex() to all edges and set them to new positions
  // 2. insert_edge() between 2 randomly selected neighboring new inserted 
  //    vertices
  // 3. insert_vertex() to the new inserted edge and set them to new positions
  // 4. insert_edge() between all other new inserted vetices of step 1 and
  //    the new inserted vertex of step 3
  // Step 1.
  eitr = p.edges_begin();
  for (int i = 0; i < num_edge; i++, ++eitr) {
    Vertex_handle vh = insert_vertex(p, eitr);
    vh->point() = edge_point_buffer[i];
  }
  Facet_iterator fitr = p.facets_begin();
  for (int i = 0; i < num_facet; i++, ++fitr) {
    // Step 2.
    Halfedge_around_facet_circulator hcir_begin = fitr->facet_begin();
    Halfedge_around_facet_circulator hcir = hcir_begin;
    
    // After linsub, the facet valence = 6
    ASSERTION_MSG(circulator_size(hcir)==6, "(ERROR) Non-triangle facet!");
    
    Halfedge_handle e1 = ++hcir;
    ++hcir;
    Halfedge_handle e2 = ++hcir;
    ++hcir;
    Halfedge_handle e3 = ++hcir;
    e2 = insert_edge(p, e1, e2);
    e3 = insert_edge(p, e2, e3);
    insert_edge(p, e3, e1);
  }

  // Update the geometry data of the newly inserted vertices by the 
  // vertices buffer
  vitr = p.vertices_begin();
  for (int i = 0; i < num_vertex; i++, ++vitr)
    vitr->point() = vertex_point_buffer[i];

  delete []vertex_point_buffer;
}


// ======================================================================
///
template <class _P> template <template <typename> class RULE>
void Polyhedron_subdivision<_P>::qt_quadralize_1step(_P& p, RULE<_P> rule) {
  p.normalize_border();

  // Build a new vertices buffer has the following structure
  //
  // 0 1 ... e_begin ... f_begin ... (end_of_buffer)
  // 0 ... e_begin-1       : store the positions of the vertex-vertices
  // e_begin ... f_begin-1 : store the positions of the edge-vertices
  // f_begin ... (end)     : store the positions of the face-vertices
  // The index of the vertices buffer should 1-1 map to the distance
  // of the corresponding iterator to the begin of the iterator.
  int num_vertex = p.size_of_vertices();
  int num_edge = p.size_of_halfedges()/2;
  int num_facet = p.size_of_facets();

  int* facet_deg = new int[num_facet];

  // If Polyhedron is using vector, we need to reserve the memory to prevent 
  // the CGAL_assertion.
  // This function for polyhedron using list is VOID.
  p.reserve(num_vertex+num_edge+num_facet, 4*2*num_edge, 4*num_edge/2);

  Point* vertex_point_buffer = new Point[num_vertex + num_edge + num_facet];
  Point* edge_point_buffer = vertex_point_buffer + num_vertex;
  Point* face_point_buffer = edge_point_buffer + num_edge;

  std::vector<bool> v_onborder(num_vertex);

  Facet_iterator fitr = p.facets_begin();
  for (int i = 0; i < num_facet; i++, ++fitr) {
    facet_deg[i] = CGAL::circulator_size(fitr->facet_begin());
    if (facet_deg[i] > 3) // Q rule
      rule.face_point_rule(fitr, face_point_buffer[i]);
  }

  int sb = p.size_of_border_edges();

  Edge_iterator eitr = p.edges_begin();
  for (int i = 0; i < num_edge-sb; i++, ++eitr)
    rule.edge_point_rule(eitr, edge_point_buffer[i]);
  for (int i = num_edge-sb; i < num_edge; i++, ++eitr) {
    int v = std::distance(p.vertices_begin(), eitr->vertex());
    v_onborder[v] = true;
    rule.border_point_rule(eitr, edge_point_buffer[i], vertex_point_buffer[v]);
  }

  Vertex_iterator vitr = p.vertices_begin();
  for (int i = 0; i < num_vertex; i++, ++vitr)
    if (!v_onborder[i]) rule.vertex_point_rule(vitr, vertex_point_buffer[i]);

  // Build the connectifty using insert_vertex() and insert_edge()
  // 1. insert_vertex() to all edges and set them to new positions
  // 2. insert_edge() between 2 randomly selected neighboring new inserted 
  //    vertices
  // 3. insert_vertex() to the new inserted edge and set them to new positions
  // 4. insert_edge() between all other new inserted vetices of step 1 and
  //    the new inserted vertex of step 3
  // Step 1.
  eitr = p.edges_begin();
  for (int i = 0; i < num_edge; i++, ++eitr) {
    Vertex_handle vh = insert_vertex(p, eitr);
    vh->point() = edge_point_buffer[i];
  }
  fitr = p.facets_begin();

  // TODO: the topoloy modification can be done by a template function
  //       and that gives the user a chance to create new topoly rules.
  for (int i = 0; i < num_facet; i++, ++fitr) {
    if (facet_deg[i] > 3) { // Q rules
      // Step 2.
      Halfedge_around_facet_circulator hcir_begin = fitr->facet_begin();
      Halfedge_around_facet_circulator hcir = hcir_begin;
      
      Halfedge_handle e1 = ++hcir; // e1 points to the newly inserted vertex
      ++hcir; // Skips one original vertex 
      Halfedge_handle e2 = ++hcir; // points to the next newly inserted vertex
      ++hcir; // Must move the cir before inserts the new edge !!
      Halfedge_handle newe = insert_edge(p, e1, e2);
      
      // Step 3.
      Halfedge_handle newv = insert_vertex_return_edge(p, newe);
      newv = newv->opposite()->prev(); // change newv to the larger face and 
                                       // still points to the newly inserted 
                                       // vertex
      // Update the geometry data of the newly inserted face-vertices
      newv->vertex()->point() = face_point_buffer[i];
   
      // Step 4.
      while (hcir != hcir_begin) {
	e1 = ++hcir;
	++hcir; // Must move the cir before inserts the new edge !!
	insert_edge(p, e1, newv); 
      }
    } else { // T rules
      // Step 2.
      Halfedge_around_facet_circulator hcir_begin = fitr->facet_begin();
      Halfedge_around_facet_circulator hcir = hcir_begin;
      
      Halfedge_handle e1 = ++hcir;
      ++hcir;
      Halfedge_handle e2 = ++hcir;
      ++hcir;
      Halfedge_handle e3 = ++hcir;
      e2 = insert_edge(p, e1, e2);
      e3 = insert_edge(p, e2, e3);
      insert_edge(p, e3, e1);
    }
  }
  
  // Update the geometry data of the newly inserted vertices by the 
  // vertices buffer
  vitr = p.vertices_begin();
  for (int i = 0; i < num_vertex; i++, ++vitr) 
    vitr->point() = vertex_point_buffer[i];

  delete []facet_deg;
  delete []vertex_point_buffer;
}

// ======================================================================
///
template <class _P> template <template <typename> class RULE>
void Polyhedron_subdivision<_P>::dualize_1step(_P& p, RULE<_P> rule) {
  int num_v = p.size_of_vertices();
  int num_e = p.size_of_halfedges()/2;
  int num_f = p.size_of_facets();
  int num_facet = num_v + num_e + num_f;
  
  // init the buffer for the next level
  Point* point_buffer = new Point[num_e*2];
  int** facet_buffer = new int*[num_facet];
  for (int i = 0; i < num_facet; ++i) facet_buffer[i] = NULL;

  // build the point_buffer
  Halfedge_iterator he_itr = p.halfedges_begin(); 
  for (int i = 0; i < num_e*2; ++i, ++he_itr) {
    Halfedge_around_facet_circulator cir = he_itr->facet_begin();
    rule.point_rule(cir, point_buffer[i]);
  }

  // build the facet_buffer
  he_itr = p.halfedges_begin(); 
  Facet_iterator fitr = p.facets_begin();  
  for (int i = 0; i < num_f; ++i, ++fitr) {
    Halfedge_around_facet_circulator  cir = fitr->facet_begin();
    int n =  CGAL::circulator_size(cir); 
    facet_buffer[i] = new int[n+1];
    facet_buffer[i][0] = n;
    for (int j = 1; j < n+1; ++j, ++cir)
      facet_buffer[i][j] = 
	std::distance(he_itr, Halfedge_handle(cir.operator->())); 
  }
  Halfedge_iterator eitr = p.halfedges_begin();
  for (int i = num_f; i < num_f+num_e; ++i, ++eitr) {
    facet_buffer[i] = new int[4+1];
    facet_buffer[i][0] = 4;
    facet_buffer[i][1] = (i-num_f)*2;
    facet_buffer[i][2] = std::distance(he_itr, eitr->prev());    
    ++eitr;
    facet_buffer[i][3] = (i-num_f)*2+1; 
    facet_buffer[i][4] = std::distance(he_itr, eitr->prev());    
  }
  Vertex_iterator vitr = p.vertices_begin();
  for (int i = num_f+num_e; i < num_f+num_e+num_v; ++i, ++vitr) {
    Halfedge_around_vertex_circulator  cir = vitr->vertex_begin();
    int n =  CGAL::circulator_size(cir); 
    facet_buffer[i] = new int[n+1];
    facet_buffer[i][0] = n;

    for (int j = 1; j < n+1; ++j, --cir)
      facet_buffer[i][j] = 
	std::distance(he_itr, Halfedge_handle(cir.operator->())); 
  }
  
  p.clear();
  Polyhedron_memory_builder<Polyhedron> pb(num_e*2, point_buffer, 
					   num_f+num_e+num_v, facet_buffer);
  p.delegate(pb);
  
  // release the buffer of the new level
  for (int i = 0; i < num_facet; ++i) delete[] facet_buffer[i];
  delete[] facet_buffer;
  delete[] point_buffer;
}



// ======================================================================
///
template <class _P> template <template <typename> class RULE>
void Polyhedron_subdivision<_P>::sqrt3refine_1step(_P& p, RULE<_P> rule) {
  //
  p.normalize_border();

  //
  int num_v = p.size_of_vertices();
  int num_e = p.size_of_halfedges()/2;
  int num_f = p.size_of_facets();

  p.reserve(num_v+num_f, (num_e+3*num_f)*2, 3*num_f);
 
  // prepare the smoothed center points
  Point* cpt = new Point[num_f]; 
  Facet_iterator fitr = p.facets_begin();
  for (int i = 0; i < num_f; ++i, ++fitr) {
    ASSERTION_MSG(circulator_size(fitr->facet_begin())==3, "(ERROR) Non-triangle facet!");
    rule.face_point_rule(fitr, cpt[i]);
  }

  // smooth the vertex points
  Vertex_iterator vitr = p.vertices_begin();
  for (int i = 0; i < num_v; ++i, ++vitr)
    rule.vertex_point_rule(vitr, vitr->point());
  
  // insert the facet points
  fitr = p.facets_begin();
  for (int i = 0; i < num_f; ++i, ++fitr) {
    Halfedge_handle center = p.create_center_vertex(fitr->halfedge());
    center->vertex()->point() = cpt[i];
  }
  
  delete []cpt;

  // flip the old edges except the border edges
  Edge_iterator eitr = p.edges_begin();
  for (int i = 0; i < num_e; ++i) {
    Halfedge_handle e = eitr;
    ++eitr; // move to next edge before flip since flip destroys current edge
    if (!e->is_border_edge()) {
      Halfedge_handle h = p.join_facet(e);
      p.split_facet(h->prev(), h->next());
    }
  }

  // TODO: border ...
  

  CGAL_postcondition(p.is_valid());
}

#endif //_POLYHEDRON_SUBDIVISION_H_01292002
