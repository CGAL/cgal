#ifndef CGAL_AABB_SEGMENT_2_PRIMITIVE_H
#define CGAL_AABB_SEGMENT_2_PRIMITIVE_H

namespace CGAL {

template <class GeomTraits, class Iterator_, class ContainerType>
class AABB_segment_2_primitive {

public:

    typedef typename GeomTraits::Point_2 Point;
    typedef typename GeomTraits::Segment_2 Datum;
    typedef ContainerType Container;
    typedef Iterator_ Id;

private:

    Id m_it;
    Datum m_datum;

public:

    AABB_segment_2_primitive() {}

    AABB_segment_2_primitive(Id it) : m_it(it) {
        m_datum = *it;
    }

    AABB_segment_2_primitive(const AABB_segment_2_primitive &primitive) {
        m_it = primitive.id();
        m_datum = primitive.datum();
    }

public:

    Id &id() {
        return m_it;
    }

    const Id &id() const {
        return m_it;
    }

    Datum &datum() {
        return m_datum;
    }

    const Datum &datum() const {
        return m_datum;
    }

    // Return a point on the primitive
    Point reference_point() const {
        return m_datum.source();
    }
};

} // namespace CGAL

#endif
