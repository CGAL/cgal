namespace CGAL {
/*!
\ingroup PkgMesh2Enum

The enum `Mesh_optimization_return_code` is the output of the global
mesh optimization functions. This output corresponds to mesh
optimization process termination reasons.  Note that depending on what
parameters have been set to the optimizer, each return value may
represent a failure or a success.

\sa `CGAL::lloyd_optimize_mesh_2`

*/
enum Mesh_optimization_return_code {
BOUND_REACHED = 0, ///< The given lower bound on mesh quality is reached.
TIME_LIMIT_REACHED, ///< The given time limit is reached.
CANT_IMPROVE_ANYMORE, ///< Mesh could not be improved anymore.
CONVERGENCE_REACHED, ///< The given convergence bound is reached.
MAX_ITERATION_NUMBER_REACHED, ///< The given maximum iteration number is reached.
ALL_VERTICES_FROZEN ///< All vertices have been frozen.
};

}
