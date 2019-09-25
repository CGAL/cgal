namespace CGAL{
//! \ingroup IOstreamFunctions
/// Reads a file with `.stl` format.
///
/// \tparam PointRange must be a model of the concept `RandomAccessContainer` or a %CGAL point type
/// \tparam TriangleRange must be a model of the concept `RandomAccessContainer`
/// with a `value_type` being a model of the concept `RandomAccessContainer`, with
/// a `value_type` begin `std::size_t`.
///
/// \param input the input stream
/// \param points a container that will contain the points used in the .stl file
/// \param facets a container that will contain the triangles used in the .stl file
/// \param verbose whether to enable or not a sanity log
///
/// \returns `true` if the reading process went well, `false` otherwise
///
/// \warning `points` and `facets` are not cleared: new points and triangles are added to the back
///          of the containers.
///
/// Although the STL file format uses triangles, it is convenient to be able to use vectors
/// and other models of the `SequenceContainer` (instead of arrays) for the face type,
/// to avoid having to convert the to apply polygon soup reparation algorithms.
template <class PointRange, class TriangleRange>
bool read_STL(std::istream& input,
              PointRange& points,
              TriangleRange& facets,
              bool verbose = false);
}
