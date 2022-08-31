#ifndef CGAL_DOMAIN_CELL_TYPE
#define CGAL_DOMAIN_CELL_TYPE

#include <limits>

namespace CGAL {
namespace Isosurfacing {

typedef std::size_t Cell_type;

static constexpr Cell_type ANY_CELL = std::numeric_limits<Cell_type>::max();

static constexpr Cell_type POLYHERDAL_CELL = ((Cell_type)1) << 0;
static constexpr Cell_type TETRAHEDRAL_CELL	= ((Cell_type)1) << 1;
static constexpr Cell_type CUBICAL_CELL = ((Cell_type)1) << 2;

}
}

#endif  // CGAL_DOMAIN_CELL_TYPE