/**
 * @file   data/3d/skel/ConstOffsetEvent.h
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#ifndef DATA_3D_SKEL_CONSTOFFSETEVENT_H
#define DATA_3D_SKEL_CONSTOFFSETEVENT_H

#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class ConstOffsetEvent : public AbstractEvent {
public:
    virtual ~ConstOffsetEvent();
    static ConstOffsetEventSPtr create();
    static ConstOffsetEventSPtr create(double offset);

    double getOffset() const;
    void setOffset(double offset);
protected:
    ConstOffsetEvent();
    ConstOffsetEvent(double offset);
    double offset_;
};

} } }

#endif /* DATA_3D_SKEL_CONSTOFFSETEVENT_H */
