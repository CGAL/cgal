namespace CGAL{
  /**
 * @brief outputs mesh to avizo format
 * @param os the stream
 * @param c3t3 the mesh
 * \see \ref IOStreamAvizo
 */
template <class C3T3>
void
output_to_avizo(std::ostream& os,
                 const C3T3& c3t3);
}
