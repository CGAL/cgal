#ifndef CGAL_TRIANGULATE_HOLE_H
#define CGAL_TRIANGULATE_HOLE_H

#include <CGAL/internal/Hole_filling/Triangulate_hole_Polyhedron_3.h>

namespace CGAL {

  /*!
  \ingroup PkgPolygonMeshProcessing
  Function triangulating a hole in a polygon mesh.
  The hole should contain no non-manifold vertex. Generated patch is guaranteed to not to break edge manifoldness and contain no degenerate triangle.
  If no possible patch is found, @a pmesh is not altered in any way, and no face descriptor is put into @a out.

  @tparam PolygonMesh must be a model of `MutableFaceGraph`
  @tparam OutputIterator iterator holding `boost::graph_traits<PolygonMesh>::face_descriptor` for patch faces.

  @param pmesh polygon mesh containing the hole
  @param border_halfedge a border halfedge incident to the hole
  @param out iterator over patch faces.
  @param use_delaunay_triangulation

  @return @a out

  \todo handle islands
  */
  template<class PolygonMesh, class OutputIterator>
  OutputIterator
    triangulate_hole(PolygonMesh& pmesh,
    typename boost::graph_traits<PolygonMesh>::halfedge_descriptor border_halfedge,
    OutputIterator out,
    bool use_delaunay_triangulation = false)
  {
    CGAL_precondition(face(border_halfedge, pmesh) == boost::graph_traits<PolygonMesh>::null_face());
    return internal::triangulate_hole_Polyhedron
      (pmesh, border_halfedge, out, use_delaunay_triangulation).first;
  }

} //end namespace CGAL

#endif //CGAL_TRIANGULATE_HOLE_H
