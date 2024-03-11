/**
 * @file   data/2d/skel/ConstOffsetEvent.h
 * @author Gernot Walzl
 * @date   2012-04-06
 */

#ifndef DATA_2D_SKEL_CONSTOFFSETEVENT_H
#define DATA_2D_SKEL_CONSTOFFSETEVENT_H

#include "debug.h"
#include "data/2d/skel/ptrs.h"
#include "data/2d/skel/AbstractEvent.h"

namespace data { namespace _2d { namespace skel {

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

#endif /* DATA_2D_SKEL_CONSTOFFSETEVENT_H */

