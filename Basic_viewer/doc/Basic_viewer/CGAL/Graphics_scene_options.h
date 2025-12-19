
namespace CGAL {

/*!
\ingroup PkgBasicViewerClasses

The class `Graphics_scene_options` is used to tune the way that the cells of a given data structure of \cgal are considered.
The different `std::function` can be modified to change for example the behavior of the drawing.
`VolumeDescriptor` can be `void` for data structures that do not represent volumes.

This class is a model of `GraphicsSceneOptions` when `VolumeDescriptor` is `void`, or a model of `GraphicsSceneOptionsWithVolumes` otherwise (`VolumeDescriptor` non `void`).

\tparam DS a data structure of \cgal.
\tparam VertexDescriptor a descriptor of vertices of `DS`.
\tparam EdgeDescriptor a descriptor of edges of `DS`.
\tparam FaceDescriptor a descriptor of faces of `DS`.
\tparam VolumeDescriptor a descriptor of volumes of `DS`. `void` by default.

\cgalModels{GraphicsSceneOptions or GraphicsSceneOptionsWithVolumes}
*/

template <typename DS,
          typename VertexDescriptor,
          typename EdgeDescriptor,
          typename FaceDescriptor,
          typename VolumeDescriptor=void>
struct Graphics_scene_options
{
public:
  typedef VertexDescriptor vertex_descriptor;
  typedef EdgeDescriptor edge_descriptor;
  typedef FaceDescriptor face_descriptor;
  typedef VolumeDescriptor volume_descriptor;

};

} // End namespace CGAL
