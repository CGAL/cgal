namespace CGAL {
  
/*!
\ingroup PkgDrawFaceGraphWithPaths

Open a new window and draw `amesh`,  either a 2D linear cell complex or a model of the FaceGraph concept, plus the paths lying on this mesh given in `paths`.
The function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam Mesh either a 2D linear cell complex or a model of the FaceGraph concept.
\param amesh the mesh to draw.
\param apaths the paths to draw, which should lie on `amesh`.
*/
template<class Mesh>
void draw(const Mesh& amesh,
          const std::vector<Path_on_surface<Mesh> >& apaths);

} /* namespace CGAL */

