namespace CGAL {

//! \ingroup IOstreamFunctions
///reads a file in .obj format.
///\tparam Points_3 a RandomAccessContainer of Point_3,
///\tparam Faces a RandomAccessContainer of RandomAccessContainer of std::size_t
/// \see IOStreamOBJ
template <class Points_3, class Faces>
bool
read_OBJ(std::istream& input, Points_3 &points, Faces &faces);
}
