namespace CGAL {

/*!
\ingroup PkgDrawFaceGraphWithPaths

opens a new window and draws `amesh`, either a 2D linear cell complex or a model of the FaceGraph concept, plus the paths lying on this mesh given in `apaths`. Parameters of the drawing are taken from the optional graphics scene options parameter.

A call to this function blocks the execution of the program until the drawing window is closed. This function requires `CGAL_Qt6`, and is only available if the macro `CGAL_USE_BASIC_VIEWER` is defined.
Linking with the cmake target `CGAL::CGAL_Basic_viewer` will link with `CGAL_Qt6` and add the definition `CGAL_USE_BASIC_VIEWER`.

\tparam Mesh either a 2D linear cell complex or a model of the FaceGraph concept.
\tparam GSOptions a model of `GraphicsSceneOptionsFaceGraphWithPaths` concept.

\param amesh the mesh to draw.
\param apaths the paths to draw, which should lie on `amesh`.
\param gso the graphics scene options parameter.
*/
template<class Mesh, class GSOptions>
void draw(const Mesh& amesh,
          const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >& apaths,
          const GSOptions& gso);

/*!
\ingroup PkgDrawFaceGraphWithPaths

A shortcut to `CGAL::draw(amesh, apaths, Graphics_scene_options_face_graph_with_paths{})`.
*/
template<class Mesh>
void draw(const Mesh& amesh,
          const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >& apaths);

/*!
\ingroup PkgDrawFaceGraphWithPaths

Same function than <a href="#gaef57f5480a700871bda6a50502338a66"><b>draw()</b></a> but taking the paths from a list instead from an std::vector.
*/
template<class Mesh, class GSOptions>
void draw(const Mesh& amesh,
          std::initializer_list<Path_on_surface<Mesh> > apaths,
          const GSOptions& gso);

/*!
\ingroup PkgDrawFaceGraphWithPaths

A shortcut to `CGAL::draw(amesh, apaths, Graphics_scene_options_face_graph_with_paths{})`.
*/
template<class Mesh>
void draw(const Mesh& amesh,
          std::initializer_list<Path_on_surface<Mesh> > apaths);

/*!
\ingroup PkgDrawFaceGraphWithPaths

adds the vertices, edges and faces of `amesh`, either a 2D linear cell complex or a model of the FaceGraph concept, plus the paths lying on this mesh given in `apaths`, into the given graphic scene `gs`. Parameters of the cells are taken from the optional graphics scene options parameter `gso`. Note that `gs` is not cleared before being filled (to enable to draw several data structures in the same basic viewer).

\tparam Mesh either a 2D linear cell complex or a model of the FaceGraph concept.
\tparam GSOptions a model of `GraphicsSceneOptionsFaceGraphWithPaths` concept.

\param amesh the mesh to draw.
\param gs the graphic scene to fill.
\param gso the graphics scene options parameter.
\param apaths the paths to draw, which should lie on `amesh`.
*/
template <class Mesh, class GSOptions>
void add_to_graphics_scene(const Mesh& amesh,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso,
                           const std::vector<Surface_mesh_topology::Path_on_surface<Mesh>>& apaths);

/*!
\ingroup PkgDrawFaceGraphWithPaths

A shortcut to `CGAL::add_to_graphics_scene(amesh, gs, Graphics_scene_options_face_graph_with_paths{}, apaths)`.
*/
template <class Mesh>
void add_to_graphics_scene(const Mesh& amesh,
                           CGAL::Graphics_scene& gs,
                           const std::vector<Surface_mesh_topology::Path_on_surface<Mesh>>& apaths);

/*!
\ingroup PkgDrawFaceGraphWithPaths

Same function than <a href="#gaa011d277d42130437c8665030ea12f97"><b>add_to_graphics_scene()</b></a> but taking the paths from a list instead from an std::vector.
*/
template <class Mesh, class GSOptions>
void add_to_graphics_scene(const Mesh& amesh,
                           CGAL::Graphics_scene& gs,
                           const GSOptions& gso,
                           std::initializer_list<Path_on_surface<Mesh> > apaths);

/*!
\ingroup PkgDrawFaceGraphWithPaths

A shortcut to `CGAL::add_to_graphics_scene(amesh, gs, Graphics_scene_options_face_graph_with_paths{}, apaths)`.
*/
template <class Mesh>
void add_to_graphics_scene(const Mesh& amesh,
                           CGAL::Graphics_scene& gs,
                           std::initializer_list<Path_on_surface<Mesh> > apaths);

} /* namespace CGAL */

