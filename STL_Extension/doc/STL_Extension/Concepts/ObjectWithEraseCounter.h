
/*!
\ingroup PkgSTLExtensionConcepts
\cgalConcept

The concept `ObjectWithEraseCounter` describes the functionalities
an object must provide so that its erase counter is updated by
a `CGAL::Compact_container` or a `CGAL::Concurrent_compact_container`.

\sa `CGAL::Compact_container`
\sa `CGAL::Concurrent_compact_container`
*/
class ObjectWithEraseCounter {
public:
/// \name Operations
/// @{
  /// Gets the value of the erase counter
  unsigned int erase_counter() const;
  /// Sets the value of the erase counter
  void set_erase_counter(unsigned int c);
  /// Increment the value of the erase counter
  void increment_erase_counter();
/// @}
};