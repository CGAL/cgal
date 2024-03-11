/**
 * @file   data/3d/skel/SphericalDblReturnEvent.h
 * @author Gernot Walzl
 * @date   2013-01-28
 */

#ifndef DATA_3D_SKEL_SPHERICALDBLRETURNEVENT_H
#define DATA_3D_SKEL_SPHERICALDBLRETURNEVENT_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/SphericalAbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SphericalDblReturnEvent : public SphericalAbstractEvent {
public:
    virtual ~SphericalDblReturnEvent();

    static SphericalDblReturnEventSPtr create();

    double getOffset() const;
    void setOffset(double offset);
    CircularVertexSPtr getVertex1() const;
    void setVertex1(CircularVertexSPtr vertex_1);
    CircularVertexSPtr getVertex2() const;
    void setVertex2(CircularVertexSPtr vertex_2);

    void setHighlight(bool highlight);

protected:
    SphericalDblReturnEvent();
    double offset_;
    CircularVertexSPtr vertex_1_;
    CircularVertexSPtr vertex_2_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALDBLRETURNEVENT_H */

