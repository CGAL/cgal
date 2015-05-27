
/*!
\ingroup PkgStlExtensionConcepts
\cgalConcept

A type `Key` is a model of the concept `Hashable` if the
specializations `boost::hash<Key>` and `std::hash<Key>` exist. 

\cgalHasModel All handles and indices of \cgal data structures.


\sa `CGAL::Unique_hash_map<Key,Data,Hash>` 
\sa `std::unordered_set`
\sa `std::unordered_map`
\sa `boost::unordered_set`
\sa `boost::unordered_map`

*/

class Hashable {

};

