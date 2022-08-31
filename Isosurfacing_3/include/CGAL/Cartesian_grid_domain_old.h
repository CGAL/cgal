#ifndef CGAL_CARTESIAN_GRID_DOMAIN_OLD_H
#define CGAL_CARTESIAN_GRID_DOMAIN_OLD_H

#include <CGAL/Cartesian_grid_3.h>

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits>
class Cartesian_grid_domain_old {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point_3;

public:
    Cartesian_grid_domain_old(const Cartesian_grid_3<Geom_traits>& grid) : grid(&grid) {}

    std::size_t size_x() const {
        return grid->xdim();
    }
    std::size_t size_y() const {
        return grid->ydim();
    }
    std::size_t size_z() const {
        return grid->zdim();
    }

    Point_3 position(const std::size_t x, const std::size_t y, const std::size_t z) const {
        const FT vx = grid->voxel_x();
        const FT vy = grid->voxel_y();
        const FT vz = grid->voxel_z();

        return Point_3(x * vx + grid->offset_x(), y * vy + grid->offset_y(), z * vz + grid->offset_z());
    }

    FT value(const std::size_t x, const std::size_t y, const std::size_t z) const {
        return grid->value(x, y, z);
    }

private:
    const Cartesian_grid_3<Geom_traits>* grid;
};

}  // namespace Isosurfacing
}  // namespace CGAL

#endif  // CGAL_CARTESIAN_GRID_DOMAIN_OLD_H
