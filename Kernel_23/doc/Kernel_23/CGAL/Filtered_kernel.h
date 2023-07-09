
namespace CGAL {

/*!
\ingroup kernel_classes

\brief `Filtered_kernel_adaptor` is a kernel that uses a filtering technique to obtain a kernel with exact and efficient predicate functors.

\details
The geometric constructions are exactly those
of the kernel `CK`, which means that they are not necessarily exact.
More details about the filtering technique can be found in  \cgalCite{cgal:bbp-iayed-01}.

In contrast to `Filtered_kernel`,
the global functions are those of `CK`.

\cgalModels{Kernel}

\cgalHeading{Example}

The following example shows how to produce a kernel whose geometric
objects and constructions are those of `Simple_cartesian<double>`
The predicate functors of the kernel are exact, the global functions
are not.

\code
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

typedef CGAL::Simple_cartesian<double> CK;
typedef CGAL::Filtered_kernel_adaptor<CK> K;

typedef K::Point_2 p(0,0), q(1,1), r(1,5);

CGAL::orientation(p,q,r); // not exact

typedef K::Orientation_2 orientation;
orientation(p,q,r); // exact
\endcode

*/
template< typename CK >
struct Filtered_kernel_adaptor {

}; /* end Filtered_kernel_adaptor */
} /* end namespace CGAL */

namespace CGAL {

/*!
\ingroup kernel_classes

`Filtered_kernel` is a kernel that uses a filtering technique based
on interval arithmetic form to achieve exact and efficient predicates.
It is based on \cgalCite{cgal:bbp-iayed-01}. In addition, a few selected important
predicates are implemented using the formally proved, semi-static, filtering
techniques from \cgalCite{cgal:mp-fcafg-05}.

The geometric constructions are exactly those
of the kernel `CK`, which means that they are not necessarily exact.

\cgalHeading{Parameters}

The first parameter, `CK`, is the "Construction Kernel", namely the kernel
from which are taken the types of the geometric objects as well as the
geometric constructions.

The second parameter, `UseStaticFilters`, is a Boolean value which
activates or not an additional layer of semi-static filters. It defaults to
`true` (activated), unless the `CGAL_NO_STATIC_FILTERS` macro is
defined. This option is mostly for debugging and testing, there should be no
production use for deactivating static filters.

\cgalModels{Kernel}

\cgalHeading{Example}

The following example shows how to produce a kernel whose geometric
objects and constructions are those of `Simple_cartesian<double>`
but the predicates are exact.

\code
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

typedef CGAL::Simple_cartesian<double> CK;
typedef CGAL::Filtered_kernel<CK> K;
\endcode

\cgalHeading{Implementation}

The implementation uses `CGAL::Filtered_predicate<EP, FP, C2E, C2F>` over
each predicate of the kernel traits interface. Additionally, faster static
filters may be used for a few selected critical predicates.

*/
template< typename CK >
struct Filtered_kernel {

/// \name Types
/// @{

/*!
The type of the exact kernel.
*/
typedef EK Exact_kernel;

/*!
The type of the approximate "filtering" kernel.
*/
typedef FK Approximate_kernel;

/// @}

/// \name Constants
/// @{

/*!
A Boolean value corresponding to the
second template argument. Tells whether static filters are provided.
*/
static const bool Has_static_filters;

/// @}

}; /* end Filtered_kernel */
} /* end namespace CGAL */
