// Copyright (c) 2025 GeometryFactory (France) and Telecom Paris (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Iasonas Manolas, Jane Tournois

#ifndef CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_REMESH_IMPL_H
#define CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_REMESH_IMPL_H

#include <CGAL/license/Tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/internal/elementary_operations.h>
#include <CGAL/Tetrahedral_remeshing/internal/edge_split_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/edge_collapse_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/edge_flip_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/vertex_smooth_operation.h>
#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>

#include <type_traits>
#include <iostream>
#include <memory>

namespace CGAL {
namespace Tetrahedral_remeshing {
namespace internal {

template <typename C3t3, typename SizingFunction, typename CellSelector, typename Visitor>
class Elementary_remesher
{
  typedef typename C3t3::Triangulation Tr;
  typedef typename Tr::Geom_traits::FT FT;

  typedef typename Tr::Cell_handle Cell_handle;
  typedef typename Tr::Vertex_handle Vertex_handle;
  typedef typename Tr::Edge Edge;

#if defined CGAL_CONCURRENT_TETRAHEDRAL_REMESHING && defined CGAL_LINKED_WITH_TBB
  template <typename Operation>
  using ExecutionPolicy = std::conditional_t<
        std::is_convertible<typename Tr::Concurrency_tag, CGAL::Parallel_tag>::value,
        ElementaryOperationExecutionParallel<Operation>,
        ElementaryOperationExecutionSequential<Operation> >;
#else
  template <typename Operation>
  using ExecutionPolicy = ElementaryOperationExecutionSequential<Operation>;
#endif

public:

  typedef EdgeSplitOperation<C3t3, SizingFunction, CellSelector> EdgeSplitOp;
  typedef EdgeCollapseOperation<C3t3, SizingFunction, CellSelector,Visitor> EdgeCollapseOp;
  typedef InternalEdgeFlipOperation<C3t3, CellSelector, Visitor> InternalEdgeFlipOp;
  typedef BoundaryEdgeFlipOperation<C3t3, CellSelector, Visitor> BoundaryEdgeFlipOp;
  typedef InternalVertexSmoothOperation<C3t3, SizingFunction, CellSelector> InternalVertexSmoothOp;
  typedef SurfaceVertexSmoothOperation<C3t3, SizingFunction, CellSelector> SurfaceVertexSmoothOp;
  typedef ComplexEdgeVertexSmoothOperation<C3t3, SizingFunction, CellSelector> ComplexEdgeSmoothOp;

  using VertexSmoothingContext = VertexSmoothingContext<C3t3, SizingFunction, CellSelector>;

  std::unique_ptr<ComplexEdgeSmoothOp> m_edge_smooth_op;
  std::unique_ptr<SurfaceVertexSmoothOp> m_surface_vertices_smooth_op;
  std::unique_ptr<InternalVertexSmoothOp> m_internal_vertices_smooth_op;
  std::shared_ptr<VertexSmoothingContext> m_context;
  Elementary_remesher()
    : m_edge_smooth_op(nullptr), m_surface_vertices_smooth_op(nullptr), m_internal_vertices_smooth_op(nullptr), m_context(nullptr), m_visitor(nullptr)
  {}

  Elementary_remesher(C3t3& c3t3, const SizingFunction& sizing, const CellSelector& cell_selector, const Visitor* visitor = nullptr)
    : m_edge_smooth_op(nullptr), m_surface_vertices_smooth_op(nullptr), m_internal_vertices_smooth_op(nullptr), m_context(nullptr), m_visitor(visitor)
  {}

private:
  const Visitor* m_visitor;

public:
  static void split(C3t3& c3t3, const SizingFunction& sizing, const CellSelector& cell_selector, const bool protect_boundaries) {
    EdgeSplitOp split_op( sizing, cell_selector, protect_boundaries);
    ExecutionPolicy<EdgeSplitOp> executor;
    executor.execute(split_op, c3t3);
  }

  static void collapse(C3t3& c3t3, const SizingFunction& sizing, const CellSelector& cell_selector,const Visitor& visitor, const bool protect_boundaries) {
     EdgeCollapseOp collapse_op( sizing, cell_selector,protect_boundaries,visitor);
     ExecutionPolicy<EdgeCollapseOp> executor;
     executor.execute(collapse_op, c3t3);
  }

  static void flip(C3t3& c3t3, CellSelector& cell_selector,  Visitor& visitor, const bool protect_boundaries) {
     // Flip internal edges
     InternalEdgeFlipOp internal_flip_op(c3t3, cell_selector, protect_boundaries, visitor);
     ExecutionPolicy<InternalEdgeFlipOp> internal_executor;
     internal_executor.execute(internal_flip_op, c3t3);

     // Flip boundary edges if not protecting boundaries
     if (!protect_boundaries) {
       BoundaryEdgeFlipOp boundary_flip_op(c3t3, cell_selector, protect_boundaries, visitor);
       ExecutionPolicy<BoundaryEdgeFlipOp> boundary_executor;
       boundary_executor.execute(boundary_flip_op, c3t3);
     }
  }

  void smooth_init(C3t3& c3t3,
                   const SizingFunction& sizing,
                   const CellSelector& cell_selector,
                   const bool protect_boundaries,
                   const bool smooth_constrained_edges)
  {
      // Create shared context for vertex smoothing operations
      m_context=std::make_shared<VertexSmoothingContext>(c3t3, cell_selector, protect_boundaries,smooth_constrained_edges);
      m_edge_smooth_op = std::make_unique<ComplexEdgeSmoothOp>(c3t3, sizing, cell_selector, protect_boundaries,smooth_constrained_edges,m_context);
      m_surface_vertices_smooth_op = std::make_unique<SurfaceVertexSmoothOp>(c3t3, sizing, cell_selector, protect_boundaries,smooth_constrained_edges,m_context);
      m_internal_vertices_smooth_op =
          std::make_unique<InternalVertexSmoothOp>(c3t3, sizing, cell_selector, protect_boundaries, smooth_constrained_edges, m_context);
  }

  void smooth(
      C3t3& c3t3,
      const SizingFunction& sizing,
      const CellSelector& cell_selector,
      const bool protect_boundaries,
      const bool smooth_constrained_edges)
  {
      // Refresh context data right before use to ensure triangulation hasn't changed
      assert(m_context);
      m_context->refresh(c3t3);

      if(!protect_boundaries) {
          if(smooth_constrained_edges && m_edge_smooth_op) {
              ExecutionPolicy<ComplexEdgeSmoothOp> executor;
              executor.execute(*m_edge_smooth_op, c3t3);
          }
          #ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
          else if(!m_edge_smooth_op) {
            std::cerr << "Complex edge smoothing operation not initialized." << std::endl;
          }
          #endif

          // Smooth vertices on surface
          if(m_surface_vertices_smooth_op) {
              ExecutionPolicy<SurfaceVertexSmoothOp> executor;
              executor.execute(*m_surface_vertices_smooth_op, c3t3);
          }
          #ifdef CGAL_TETRAHEDRAL_REMESHING_VERBOSE
          else{
            std::cerr << "Surface vertex smoothing operation not initialized." << std::endl;
          }
          #endif
      }

      // Smooth internal vertices
      {
          ExecutionPolicy<InternalVertexSmoothOp> executor;
          executor.execute(*m_internal_vertices_smooth_op, c3t3);
      }
  }
  };


} // namespace internal
} // namespace Tetrahedral_remeshing
} // namespace CGAL

#endif // CGAL_TETRAHEDRAL_REMESHING_ELEMENTARY_REMESH_IMPL_H