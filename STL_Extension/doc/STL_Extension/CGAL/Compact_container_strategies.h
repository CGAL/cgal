
namespace CGAL {

/// \ingroup PkgStlExtension

/*!
\ingroup CompactContainer

The class `Compact_container_strategy_base` is a model of the concept 
`CompactContainerStrategy`. It is a basic \"do nothing\" strategy.

\cgalModels `CompactContainerStrategy`

\sa `Compact_container` 
\sa `Concurrent_compact_container`
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

\sa `Compact_container`
\sa `Concurrent_compact_container`
\sa `CC_safe_handle`
*/


// A strategy managing an internal counter
class Compact_container_strategy_with_counter
{
public:
};

/*!
\ingroup CompactContainer

The class `CC_safe_handle` is a helper class which embed a
`Compact_container` (or `Concurrent_compact_container`) iterator and 
is able to know if the pointee has been erased since the creation of the
`CC_safe_handle` instance. 
It stores the iterator and the erase counter value of the pointee when the 
`CC_safe_handle` instance was created.
It can only be used when the compact container strategy uses an erase counter.

\tparam CC_iterator is the type of the embedded `Compact_container` 
                    (or `Concurrent_compact_container`) iterator

\sa `Compact_container` 
\sa `Concurrent_compact_container` 
\sa `CompactContainerStrategy` 
*/

template <typename CC_iterator>
class CC_safe_handle
{
  typedef typename CC_iterator::Strategy Strategy;

public:
  /// \name Creation 
  /*! 
  introduces a safe handle from a `Compact_container` iterator
  */ 
  CC_safe_handle(CC_iterator handle);
  
  /// \name Access Member Functions 
  /// @{ 
  /*!
  Returns true if the pointee has been erased, i.e. if the iterator points to
  a freed object or to another object.
  */
  bool is_zombie() const;
  
  /*!
  Returns the embedded `Compact_container` iterator.
  */
  CC_iterator get_cc_iterator() const;
  /// @} 
};


/*!
\ingroup CompactContainer

The class `make_cc_safe_handle` function allows to build a `CC_safe_handle`
from a `Compact_container` (or `Concurrent_compact_container`) iterator.

\sa `CC_safe_handle`
*/

template <typename CC_iterator>
CC_safe_handle<CC_iterator> make_cc_safe_handle(CC_iterator handle);

} //namespace CGAL