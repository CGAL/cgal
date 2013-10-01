
/*!
\ingroup PkgStlExtensionConcepts
\cgalConcept

The concept `ObjectWithEraseCounter` is intended to be used by objects stored
in a `CGAL::Compact_container` or a `CGAL::Concurrent_compact_container`.

\sa `CGAL::Compact_container`
\sa `CGAL::Concurrent_compact_container`
*/
class ObjectWithEraseCounter {
public:
/// \name Operations 
/// @{
  /// Sets the value of the erase counter
  void set_erase_counter(unsigned int c);
  /// Increment the value of the erase counter
  void increment_erase_counter();
/// @}
};