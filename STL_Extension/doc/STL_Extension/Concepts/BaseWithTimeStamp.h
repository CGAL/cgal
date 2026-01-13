
/*!
\ingroup PkgSTLExtensionConcepts
\cgalConcept

The concept `BaseWithTimeStamp` describes the functionalities
an object must provide so that its timestamp is updated by
a `CGAL::Compact_container` or a `CGAL::Concurrent_compact_container`.

\cgalHasModelsBegin
\cgalHasModels{CGAL::Triangulation_simplex_base_with_time_stamp}
\cgalHasModelsEnd

\sa `CGAL::Compact_container`
\sa `CGAL::Concurrent_compact_container`
*/
class BaseWithTimeStamp {
public:
  /// Tag type to indicate that the class has a time stamp.
  /// Must be `CGAL::Tag_true`.
  using Has_timestamp = CGAL::Tag_true;

  /// @name Access Functions
  ///
  /// Member functions to read and set the time stamp.
  /// @{
  std::size_t time_stamp() const{
  }
  void set_time_stamp(const std::size_t&) {
  }
  ///@}
};
