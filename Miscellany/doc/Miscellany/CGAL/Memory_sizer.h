
namespace CGAL {

/*!
\ingroup MiscellanyRef

The class `Memory_sizer` allows to measure the memory size used by the process.
Both the virtual memory size and the resident size are available (the resident
size does not account for swapped out memory nor for the memory which is not
yet paged-in).

\cgalHeading{Implementation}

Accessing this information requires the use of non-portable code.
Currently, there is support for Linux platforms, the Microsoft and Intel
compiler on Windows, as well as Mac OS X. If a platform is not supported, the
two member functions return 0.

*/

struct Memory_sizer {

/// \name Types
/// @{

/*!

*/
typedef std::size_t size_type;

/// @}

/// \name Creation
/// @{

/*!
%Default constructor.
*/
Memory_sizer();

/// @}

/// \name Operations
/// @{

/*!
Returns the virtual memory size in bytes.
*/
size_type virtual_size() const;

/*!
Returns the resident memory size in bytes.
*/
size_type resident_size() const;

/// @}

}; /* end Memory_sizer */
} /* end namespace CGAL */
