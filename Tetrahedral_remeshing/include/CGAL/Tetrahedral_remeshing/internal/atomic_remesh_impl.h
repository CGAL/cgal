#ifndef CGAL_TETRAHEDRAL_REMESHING_ATOMIC_REMESH_IMPL_H
#define CGAL_TETRAHEDRAL_REMESHING_ATOMIC_REMESH_IMPL_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/atomic_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/edge_split_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/edge_collapse_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/edge_flip_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/vertex_smooth_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <iostream>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

template <typename C3t3, typename SizingFunction, typename CellSelector> class Atomic_remesher
{
  typedef typename C3t3::Triangulation Tr;
  typedef typename Tr::Geom_traits::FT FT;

  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Edge Edge;

  typedef EdgeSplitOperation<C3t3, SizingFunction, CellSelector> EdgeSplitOp;
  typedef EdgeCollapseOperation<C3t3, SizingFunction, CellSelector> EdgeCollapseOp;
  typedef EdgeFlipOperation<C3t3, CellSelector> EdgeFlipOp;
  typedef VertexSmoothOperation<C3t3, SizingFunction, CellSelector, SmoothingDomain::INTERNAL_VERTICES>
      InternalVertexSmoothOp;
  typedef VertexSmoothOperation<C3t3, SizingFunction, CellSelector, SmoothingDomain::SURFACE_VERTICES>
      SurfaceVertexSmoothOp;
  typedef VertexSmoothOperation<C3t3, SizingFunction, CellSelector, SmoothingDomain::COMPLEX_EDGES> ComplexEdgeSmoothOp;

  template <typename Operation> using ExecutionPolicy = AtomicOperationExecutionSequential<Operation>;

public:

  static void split(C3t3& c3t3, const SizingFunction& sizing, const CellSelector& cell_selector, const bool protect_boundaries) {
    EdgeSplitOp split_op(c3t3, sizing, cell_selector, protect_boundaries);
    ExecutionPolicy<EdgeSplitOp> executor;
    executor.execute(split_op, c3t3);
  }

  static void collapse(C3t3& c3t3, const SizingFunction& sizing, const CellSelector& cell_selector) {
     EdgeCollapseOp collapse_op(c3t3, sizing, cell_selector);
     ExecutionPolicy<EdgeCollapseOp> executor;
     executor.execute(collapse_op, c3t3);
  }

	static void flip(C3t3& c3t3, const CellSelector& cell_selector) {
	   EdgeFlipOp flip_op(c3t3, cell_selector);
	   ExecutionPolicy<EdgeFlipOp> executor;
	   executor.execute(flip_op, c3t3);
	}

  static void smooth(
      C3t3& c3t3,
      const SizingFunction& sizing,
      const CellSelector& cell_selector,
      const bool protect_boundaries,
      const bool smooth_constrained_edges)
  {
      if(!protect_boundaries) {
          // Smooth vertices on complex edges
          if(smooth_constrained_edges) {
              ComplexEdgeSmoothOp edge_smooth_op(c3t3, sizing, cell_selector, protect_boundaries, smooth_constrained_edges);
              ExecutionPolicy<ComplexEdgeSmoothOp> executor;
              executor.execute(edge_smooth_op, c3t3);
          }

          // Smooth vertices on surface
          {
              SurfaceVertexSmoothOp surface_smooth_op(c3t3, sizing, cell_selector, protect_boundaries, smooth_constrained_edges);
              ExecutionPolicy<SurfaceVertexSmoothOp> executor;
              executor.execute(surface_smooth_op, c3t3);
          }
      }

      // Smooth internal vertices
      {
          InternalVertexSmoothOp internal_smooth_op(c3t3, sizing, cell_selector, protect_boundaries, smooth_constrained_edges);
          ExecutionPolicy<InternalVertexSmoothOp> executor;
          executor.execute(internal_smooth_op, c3t3);
      }
  }
  };


} // namespace internal
} // namespace internal
} // namespace Tetrahedral_remeshing

#endif // CGAL_TETRAHEDRAL_REMESHING_ATOMIC_REMESH_IMPL_H