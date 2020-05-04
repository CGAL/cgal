namespace CGAL {

/*!
  \ingroup PkgBGLIOFct

 * \brief writes the triangle meshes contained in `gs` into the 3mf file `file_name`.
 *
 * \tparam FaceGraphRange a model of the concepts `RandomAccessContainer`
 * and `BackInsertionSequence` whose `value type` is
 * a model of the concepts `FaceListGraph` and `HalfedgeListGraph`
 * that has only triangle faces.
 *
 * \param file_name the name of the 3mf file to write.
 * \param gs a `FaceGraphRange` that contains the meshes
 *  to write. An internal property map for `CGAL::vertex_point_t`
 * must be available for each mesh.
 * \param names will contains the name of each mesh in `file_name`.
 *
 * \return `true` if the writing is successful, `false` otherwise.
 */
template<typename FaceGraphRange>
bool write_triangle_meshes_to_3mf(const std::string& file_name,
                                  const FaceGraphRange& gs,
                                  const std::vector<std::string>& names);
} // namespace CGAL
