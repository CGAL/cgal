
namespace CGAL {

/// \ingroup PkgStlExtension

/*!
\ingroup CompactContainer

The class `Compact_container_strategy_base` is a model of the concept 
`CompactContainerStrategy`. It is a basic "do nothing" strategy.

\cgalModels `CompactContainerStrategy`

\sa `CGAL::Compact_container` 
\sa `CGAL::Concurrent_compact_container` 
*/

class Compact_container_strategy_base {
public:
};


/*!
\ingroup CompactContainer

The class `Compact_container_strategy_with_counter` is a model of the concept 
`CompactContainerStrategy`. It is a strategy for elements containing an erase
counter. It can be used when the elements themselves provide three member 
functions: get_erase_counter, set_erase_counter and increment_erase_counter.

\cgalModels `CompactContainerStrategy`

\sa `CGAL::Compact_container` 
\sa `CGAL::Concurrent_compact_container` 
*/


// A strategy managing an internal counter
class Compact_container_strategy_with_counter
{
public:
};

} //namespace CGAL