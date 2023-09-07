
/*!
\ingroup PkgSTLExtensionConcepts
\cgalConcept

A type `Key` is a model of the concept `Hashable` if the
specializations `boost::hash<Key>` and `std::hash<Key>` exist.

\cgalHasModelsBegin
\cgalHasModelsBare{All handles and indices of \cgal data structures}
\cgalHasModelsBare{All handles of OpenMesh\, by including the specializations of the `graph_traits` header files provided by \cgal. They can be disables by defining the macro `CGAL_DISABLE_HASH_OPENMESH`.}
\cgalHasModelsEnd

\sa `CGAL::Unique_hash_map<Key,Mapped,Hash>`
\sa <A HREF="https://en.cppreference.com/w/cpp/container/unordered_set">`std::unordered_set`</a>
\sa <A HREF="https://en.cppreference.com/w/cpp/container/unordered_map">`std::unordered_map`</a>
\sa <A HREF="https://www.boost.org/libs/unordered/doc/html/unordered.html#unordered_set">`boost::unordered_set`</a>
\sa <A HREF="https://www.boost.org/libs/unordered/doc/html/unordered.html#unordered_map">`boost::unordered_map`</a>

*/

class Hashable {

};
