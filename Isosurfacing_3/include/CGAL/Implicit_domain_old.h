#ifndef CGAL_IMPLICIT_DOMAIN_OLD_H
#define CGAL_IMPLICIT_DOMAIN_OLD_H

#include <CGAL/Bbox_3.h>

namespace CGAL {
namespace Isosurfacing {

template <class GeomTraits, typename Function>
class Implicit_domain_old {
public:
    typedef GeomTraits Geom_traits;
    typedef typename Geom_traits::FT FT;
    typedef typename Geom_traits::Point_3 Point_3;
    typedef typename Geom_traits::Vector_3 Vector_3;

public:
    Implicit_domain_old(const Function& func, const CGAL::Bbox_3& domain, const Vector_3& resolution)
        : func(func), bbox(domain), resolution(resolution) {

        sizes[0] = domain.x_span() / resolution.x();
        sizes[1] = domain.y_span() / resolution.y();
        sizes[2] = domain.z_span() / resolution.z();
    }

    std::size_t size_x() const {
        return sizes[0];
    }
    std::size_t size_y() const {
        return sizes[1];
    }
    std::size_t size_z() const {
        return sizes[2];
    }

    Point_3 position(const std::size_t x, const std::size_t y, const std::size_t z) const {
        return Point_3(x * resolution.x() + bbox.xmin(), y * resolution.y() + bbox.ymin(),
                       z * resolution.z() + bbox.zmin());
    }

    FT value(const std::size_t x, const std::size_t y, const std::size_t z) const {
        return func(position(x, y, z));
    }

private:
    Function func;

    CGAL::Bbox_3 bbox;
    Vector_3 resolution;

    std::array<std::size_t, 3> sizes;
};


template <typename Function, class GeomTraits = typename Function::Geom_traits>
Implicit_domain_old<GeomTraits, Function> create_implicit_domain_old(const Function& func, const CGAL::Bbox_3& domain,
                                                                     const typename GeomTraits::Vector_3& resolution) {
    return Implicit_domain_old<GeomTraits, Function>(func, domain, resolution);
}

}  // namespace Isosurfacing
}  // end namespace CGAL

#endif  // CGAL_IMPLICIT_DOMAIN_OLD_H
