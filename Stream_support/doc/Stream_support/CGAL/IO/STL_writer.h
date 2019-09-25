namespace CGAL{
//! \ingroup IOstreamFunctions
/// Writes a soup into a file in the `.stl` format.
///
/// \tparam PointRange must be a model of the concept `RandomAccessContainer` or a %CGAL point type
/// \tparam TriangleRange must be a model of the concept `RandomAccessContainer`
/// with a `value_type` being a model of the concept `RandomAccessContainer`, with
/// a `value_type` begin `std::size_t`.
///
/// \param out the output stream
/// \param points a container that contains the points of the soup.
/// \param facets a container that contains the triangles of the soup
///
template <class PointRange, class TriangleRange>
std::ostream&
write_STL(const PointRange& points, const TriangleRange& facets, std::ostream& out);
}
