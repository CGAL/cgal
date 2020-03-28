
/*!
\ingroup PkgSTLExtensionConcepts
\cgalConcept

A type `Key` is a model of the concept `Hashable` if the
specializations `boost::hash<Key>` and `std::hash<Key>` exist.

\cgalHasModel All handles and indices of \cgal data structures.
\cgalHasModel All handles of OpenMesh, by including the specializations
of the `graph_traits` header files provided by \cgal.
They can be disables by defining the macro `CGAL_DISABLE_HASH_OPENMESH`.

\sa `CGAL::Unique_hash_map<Key,Mapped,Hash>`
\sa <A HREF="http://www.cplusplus.com/reference/unordered_set/unordered_set/">`std::unordered_set`</a>
\sa <A HREF="http://www.cplusplus.com/reference/unordered_set/unordered_map/">`std::unordered_map`</a>
\sa <A HREF="https://www.boost.org/doc/libs/release/doc/html/boost/unordered_set.html">`boost::unordered_set`</a>
\sa <A HREF="https://www.boost.org/doc/libs/release/doc/html/boost/unordered_map.html">`boost::unordered_map`</a>

*/

class Hashable {

};

