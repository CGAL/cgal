
/*!
\ingroup PkgStlExtensionConcepts
\cgalConcept

The concept `CompactContainerStrategy` contains static functions so that the 
`CGAL::Compact_container` (or the `CGAL::Concurrent_compact_container`) can update the
so-called erase counter stored inside the objets it contains. Thus, when an object is
erased from the container, the erase counter will be incremented. 
Note that this is meaningful only because the `CGAL::Compact_container` doesn't 
deallocate elements until the destruction or clear() of the container.
For example, this counter is used by
the 3D mesh generation engine to lazily manage the queues of bad cells:
an element in the queue is a pair containing an cell iterator and the 
erase counter value of the cell when it has been inserted. When an element
is popped from the queue, the algorithm check if the current value of
the erase counter matches the stored value. If it doesn't match, it means
the cell has been destroyed in the meantime and it ignores it. Without this
lazy management, each time a cell is destroyed, the algorithm has to look
for it in the queue and remove it. This is even more useful for the parallel 
version of the meshing process, since each thread has its own queue and looking
for a cell in all the queues would be very slow.

\cgalHasModel `CGAL::Compact_container_strategy_base` 
\cgalHasModel `CGAL::Compact_container_strategy_with_counter` 

\sa `CGAL::Compact_container` 
\sa `CGAL::Concurrent_compact_container` 

*/
class CompactContainerStrategy {
public:
/// \name Constants
/// @{
/*! 
Constant saying if this strategy uses an erase counter or not
*/ 
  static const bool Uses_erase_counter;
  
/// @} 

/// \name Static operations 
/// @{ 
  
/*! 
Return the value of the erase counter stored inside `element`
*/ 
template <typename Element>
static unsigned int get_erase_counter(const Element &element);

/*! 
Set the value of the erase counter stored inside `element` to `counter_value`
*/ 
template <typename Element>
static void set_erase_counter(Element &element, unsigned int counter_value);

/*! 
Increment the value of the erase counter stored inside `element`
*/
template <typename Element>
static void increment_erase_counter(Element &element);

/// @} 
};