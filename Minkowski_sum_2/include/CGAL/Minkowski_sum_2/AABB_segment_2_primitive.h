#ifndef CGAL_AABB_SEGMENT_2_PRIMITIVE_H
#define CGAL_AABB_SEGMENT_2_PRIMITIVE_H

namespace CGAL {

// Wraps around a Segment_2 and provides its iterator as Id
template <class GeomTraits, class Iterator_, class ContainerType>
class AABB_segment_2_primitive {

public:

    typedef Iterator_ Id;
    typedef typename GeomTraits::Segment_2 Datum;
    typedef typename GeomTraits::Point_2 Point;
    typedef ContainerType Container;

    AABB_segment_2_primitive() {}

    AABB_segment_2_primitive(Id it) : m_it(it) {
    }

    AABB_segment_2_primitive(const AABB_segment_2_primitive &primitive) {
        m_it = primitive.id();
    }

    const Id &id() const {
        return m_it;
    }

    const Datum datum() const {
        return *m_it;
    }

    // Return a point on the primitive
    Point reference_point() const {
        return datum().source();
    }

private:

    Id m_it;

};

} // namespace CGAL

#endif
