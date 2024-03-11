/**
 * @file   data/3d/SphericalPolygon.h
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#ifndef DATA_3D_SPHERICALPOLYGON_H
#define DATA_3D_SPHERICALPOLYGON_H

#include "typedefs_thread.h"
#include "data/3d/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _3d {

class SphericalPolygon : public std::enable_shared_from_this<SphericalPolygon> {
public:
    virtual ~SphericalPolygon();

    static SphericalPolygonSPtr create(Sphere3SPtr sphere);
    static SphericalPolygonSPtr create(Sphere3SPtr sphere,
            unsigned int num_edges, CircularEdgeSPtr edges[]);

    SphericalPolygonSPtr clone() const;

    Sphere3SPtr getSphere() const;

    void addVertex(CircularVertexSPtr vertex);
    bool removeVertex(CircularVertexSPtr vertex);
    void addEdge(CircularEdgeSPtr edge);
    bool removeEdge(CircularEdgeSPtr edge);

    void clear();

    SharedMutex& mutex();
    std::list<CircularVertexSPtr>& vertices();
    std::list<CircularEdgeSPtr>& edges();

    bool isConsistent() const;

    double getRadius() const;

    std::string toString() const;

protected:
    SphericalPolygon(Sphere3SPtr sphere);
    mutable SharedMutex mutex_;
    Sphere3SPtr sphere_;
    std::list<CircularVertexSPtr> vertices_;
    std::list<CircularEdgeSPtr> edges_;
};

} }

#endif /* DATA_3D_SPHERICALPOLYGON_H */
