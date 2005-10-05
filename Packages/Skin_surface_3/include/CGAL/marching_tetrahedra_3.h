#ifndef MARCHING_TETRAHEDRA_H
#define MARCHING_TETRAHEDRA_H

#define COARSE_PERTURBATION  1e-05

#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Modifier_base.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Marching_tetrahedra_observer_default_3.h>

CGAL_BEGIN_NAMESPACE 

// If TriangulationDataStructure_3 only gets a Cell_handle range
// it is not possible to derive the Vertex_handle type.
template <class Triangulation_3,
	  class HalfedgeDS_3,
	  class MarchingTetrahedraTraits_3,
	  class MarchingTetrahedraObserver_3 >
class Marching_tetrahedra_builder : public Modifier_base<HalfedgeDS_3> {
public:
//   typedef TriangulationCellIterator             T_cell_iterator;
  typedef Triangulation_3                       Triangulation;
  typedef HalfedgeDS_3                          HDS;
  typedef MarchingTetrahedraTraits_3            Traits;
  typedef MarchingTetrahedraObserver_3          Observer;

private:
//   typedef typename T_cell_iterator::value_type          T_cell;
//   typedef typename T_cell::Vertex_handle                T_vertex_handle;
  typedef typename Triangulation::Finite_cells_iterator T_cell_iterator;
  typedef typename Triangulation::Vertex_handle         T_vertex_handle;

  typedef typename HDS::Face_handle                     HDS_face_handle;
  typedef typename HDS::Vertex_handle                   HDS_vertex_handle;

  typedef std::map<T_vertex_handle,bool>                T_vertex_map;
  typedef typename T_vertex_map::iterator               T_vertex_map_it;

  // First vertex lies inside the surface, the second vertex outside
  typedef std::pair<T_vertex_handle,T_vertex_handle>    T_edge;
  typedef std::map<T_edge,int>                          T_edge_map;
  typedef typename T_edge_map::iterator                 T_edge_map_it;

  typedef Polyhedron_incremental_builder_3<HDS>         Polyh_incr_builder;
public:
  Marching_tetrahedra_builder(
//     TriangulationCellIterator first, TriangulationCellIterator last,
    const Triangulation &triang,
    const Traits &traits, Observer &observer)
    : triang(triang), //first(first), last(last),
      traits(traits), observer(observer),
      nVertices(0) {
  }

  void operator()( HDS& hds) {
    // sortedV is an array of vertex indices.
    // The first nIn vertices lie inside the surface, the last 4-nIn outside.
    int sortedV[4], cellV[4], nIn;
    T_edge edge;

    Polyh_incr_builder B( hds, true);
    B.begin_surface(0,0,0);

    for (T_cell_iterator cit = triang.finite_cells_begin();
	 cit != triang.finite_cells_end(); cit++) {
      // Compute signs on vertices and sort them:
      nIn = 0;
      for (int i=0; i<4; i++) {
        if (is_inside(cit,i)) {
          sortedV[nIn] = i; nIn++;
        } else {
          sortedV[3-i+nIn] = i;
        }
      }

      // Process edges whose vertices lie on different sides of the surface
      int edgeN=0;
      for (int i=0; i<nIn; i++) {
	for (int j=nIn; j<4; j++) {
	  cellV[edgeN++] = process_edge(B, cit,sortedV[i],sortedV[j]);
	}
      }

      // Construct triangles:
      CGAL_assertion(!B.error());
      if (nIn==1) {
	process_cell(B, cellV, (sortedV[0]%2)==1, cit);
	CGAL_assertion(!B.error());
      } else if (nIn==2) {
	bool change_orientation =
	  (((sortedV[0] == 0) && (sortedV[1]==1)) || 
	    ((sortedV[0] == 0) && (sortedV[1]==3)) || 
	    ((sortedV[0] == 1) && (sortedV[1]==2)) || 
	    ((sortedV[0] == 2) && (sortedV[1]==3)));
	
	process_cell(B, cellV, !change_orientation, cit);
 	process_cell(B, cellV+1, change_orientation, cit);
	CGAL_assertion(!B.error());
      } else if (nIn==3) {
	process_cell(B, cellV, (sortedV[3]%2) == 1, cit);
	CGAL_assertion(!B.error());
      }
      
    }

    CGAL_assertion(!B.error());
    B.end_surface();
  }


  bool is_inside(T_cell_iterator ch, int i) {
    T_vertex_map_it it = triang_vertex_signs.find(ch->vertex(i));
    
    if (it == triang_vertex_signs.end()) {
      bool side = (traits.sign(ch,i) == POSITIVE);
      triang_vertex_signs[ch->vertex(i)] = side;
      CGAL_assertion(triang_vertex_signs[ch->vertex(i)] == side);
      return side;
    } else {
      return it->second;
    }
  }
  int process_edge(Polyh_incr_builder &B,
    T_cell_iterator ch, int i, int j) {
    CGAL_assertion(is_inside(ch, i));
    CGAL_assertion(!is_inside(ch, j));

    T_edge edge = T_edge(ch->vertex(i),ch->vertex(j));
    T_edge_map_it edge_it = polyh_vert.find(edge);

    if (edge_it == polyh_vert.end()) {
      HDS_vertex_handle vh = B.add_vertex(traits.intersection(ch, i, j));
      polyh_vert[edge] = nVertices;
      nVertices ++;
      observer.after_vertex_insertion(ch, i, j, vh);
      return nVertices-1;
    } else {
      return edge_it->second;
    }
  }
  
  // Orientation is right
  void process_cell(
    Polyh_incr_builder &B,
    int *vs,
    bool change_orientation,
    T_cell_iterator ch)
  {
    CGAL_assertion((vs[0]!=vs[1]) && (vs[0]!=vs[2]) && (vs[1]!=vs[2]));
    HDS_face_handle f = B.begin_facet();
    if (change_orientation) {
      B.add_vertex_to_facet( vs[0] );
      B.add_vertex_to_facet( vs[2] );
      B.add_vertex_to_facet( vs[1] );
    } else {
      B.add_vertex_to_facet( vs[0] );
      B.add_vertex_to_facet( vs[1] );
      B.add_vertex_to_facet( vs[2] );
    }
    B.end_facet();
    observer.after_facet_insertion(ch, f);
  }

private:
//   T_cell_iterator first, last;
  const Triangulation &triang;
  const Traits &traits;
  Observer &observer;
  T_edge_map polyh_vert;
  T_vertex_map triang_vertex_signs;
  int nVertices;
};

// template <class Triangulation_3,
// 	  class Polyhedron_3,
// 	  class MarchingTetrahedraTraits_3,
// 	  class MarchingTetrahedraObserver_3>
// void marching_tetrahedra_3(
//   const Triangulation_3 &triangulation,
//   Polyhedron_3 &polyhedron,
//   const MarchingTetrahedraTraits_3 &traits,
//   MarchingTetrahedraObserver_3 &observer) {

//   typedef typename Triangulation_3::Cell_iterator             Cell_iterator;
//   typedef typename Polyhedron_3::HalfedgeDS                   HDS;
//   typedef MarchingTetrahedraTraits_3                          Traits;
//   typedef MarchingTetrahedraObserver_3                        Observer;
//   typedef Marching_tetrahedra_builder<Cell_iterator,HDS, Traits, Observer>
//                                                               Builder;
  
//   Builder builder(
//     triangulation.finite_cells_begin(),
//     triangulation.finite_cells_end(),
//     traits, observer);
//   polyhedron.delegate(builder);
// }

template <class Triangulation_3,
	  class Polyhedron_3,
	  class MarchingTetrahedraTraits_3>
void marching_tetrahedra_3(
  const Triangulation_3 &triangulation,
  Polyhedron_3 &polyhedron,
  const MarchingTetrahedraTraits_3 &traits) {
  
  typedef typename Polyhedron_3::HalfedgeDS                    HDS;
  typedef MarchingTetrahedraTraits_3                           Traits;
  typedef Marching_tetrahedra_observer_default_3<Triangulation_3, Polyhedron_3>
    Observer; 
  typedef Marching_tetrahedra_builder<Triangulation_3,HDS,Traits,Observer>
    Builder;

  Observer observer;
  Builder builder(
    triangulation,
    traits, observer);
  polyhedron.delegate(builder);
}

// template <class InputIterator,
// 	  class Polyhedron_3,
// 	  class MarchingTetrahedraTraits_3,
// 	  class MarchingTetrahedraObserver_3>
// void marching_tetrahedra_3(
//   InputIterator first, InputIterator last,
//   Polyhedron_3 &polyhedron,
//   const MarchingTetrahedraTraits_3 &traits,
//   MarchingTetrahedraObserver_3 &observer) {
  
//   typedef typename Polyhedron_3::HalfedgeDS                   HDS;
//   typedef MarchingTetrahedraTraits_3                          Traits;
//   typedef MarchingTetrahedraObserver_3                        Observer;
//   typedef Marching_tetrahedra_builder<Triangulation_3,HDS, Traits, Observer>
//                                                               Builder;
  
//   Builder builder(first, last, traits, observer);
//   polyhedron.delegate(builder);
// }

// template <class Triangulation_3,
// 	  class Polyhedron_3,
// 	  class MarchingTetrahedraTraits_3>
// void marching_tetrahedra_3(
//   const Triangulation_3 &triangulation,
//   Polyhedron_3 &polyhedron,
//   const MarchingTetrahedraTraits_3 &traits) {
  
//   typedef typename Polyhedron_3::HalfedgeDS                    HDS;
//   typedef typename Triangulation_3::Cell_handle                Cell_handle;
//   typedef MarchingTetrahedraTraits_3                           Traits;
//   typedef Marching_tetrahedra_observer_default_3<Cell_handle, Polyhedron_3>
//                                                                Observer; 
//   typedef Marching_tetrahedra_builder<Triangulation_3,HDS,Traits,Observer>
//                                                                Builder;

//   Observer observer;
//   Builder builder(triangulation, traits, observer);
//   polyhedron.delegate(builder);
// }

CGAL_END_NAMESPACE 

#endif // MARCHING_TETRAHEDRA_H
