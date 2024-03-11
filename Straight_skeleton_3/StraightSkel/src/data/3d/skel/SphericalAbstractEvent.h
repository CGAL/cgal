/**
 * @file   data/3d/skel/SphericalAbstractEvent.h
 * @author Gernot Walzl
 * @date   2012-11-28
 */

#ifndef DATA_3D_SKEL_SPHERICALABSTRACTEVENT_H
#define DATA_3D_SKEL_SPHERICALABSTRACTEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include <list>
#include <string>

namespace data { namespace _3d { namespace skel {

class SphericalAbstractEvent {
public:
    virtual ~SphericalAbstractEvent();

    SphericalPolygonSPtr getPolygonResult() const;
    void setPolygonResult(SphericalPolygonSPtr polygon);
    SphericalSkeletonSPtr getSkel() const;
    void setSkel(SphericalSkeletonSPtr skel);
    std::list<SphericalAbstractEventSPtr>::iterator getListIt() const;
    void setListIt(std::list<SphericalAbstractEventSPtr>::iterator list_it);

    virtual void setHighlight(bool highlight);
    virtual double getOffset() const = 0;  // abstract

    static const int CONST_OFFSET_EVENT = 1;
    static const int EDGE_EVENT = 2;
    static const int SPLIT_EVENT = 3;
    static const int TRIANGLE_EVENT = 4;
    static const int DBL_EDGE_EVENT = 5;
    static const int LEAVE_EVENT = 6;
    static const int RETURN_EVENT = 7;
    static const int DBL_LEAVE_EVENT = 8;
    static const int DBL_RETURN_EVENT = 9;
    static const int VERTEX_EVENT = 10;
    static const int EDGE_MERGE_EVENT = 11;
    static const int INVERSION_EVENT = 12;

    virtual int getType() const;

    virtual std::string toString() const;

protected:
    SphericalAbstractEvent();

    SphericalPolygonSPtr polygon_result_;
    SphericalSkeletonWPtr skel_;
    std::list<SphericalAbstractEventSPtr>::iterator list_it_;
    int type_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALABSTRACTEVENT_H */
