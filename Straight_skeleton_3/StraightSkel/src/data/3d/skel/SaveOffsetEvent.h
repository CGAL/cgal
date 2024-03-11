/**
 * @file   data/3d/skel/SaveOffsetEvent.h
 * @author Gernot Walzl
 * @date   2013-12-28
 */

#ifndef DATA_3D_SKEL_SAVEOFFSETEVENT_H
#define DATA_3D_SKEL_SAVEOFFSETEVENT_H

#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/AbstractEvent.h"

namespace data { namespace _3d { namespace skel {

class SaveOffsetEvent : public AbstractEvent {
public:
    virtual ~SaveOffsetEvent();
    static SaveOffsetEventSPtr create();
    static SaveOffsetEventSPtr create(double offset);

    double getOffset() const;
    void setOffset(double offset);
protected:
    SaveOffsetEvent();
    SaveOffsetEvent(double offset);
    double offset_;
};

} } }

#endif /* DATA_3D_SKEL_SAVEOFFSETEVENT_H */
