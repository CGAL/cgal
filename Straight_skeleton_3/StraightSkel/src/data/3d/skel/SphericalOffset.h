/**
 * @file   data/3d/skel/SphericalOffset.h
 * @author Gernot Walzl
 * @date   2012-12-21
 */

#ifndef DATA_3D_SKEL_SPHERICALOFFSET_H
#define DATA_3D_SKEL_SPHERICALOFFSET_H

#include "data/3d/skel/ptrs.h"

namespace data { namespace _3d { namespace skel {

class SphericalOffset {
public:
    virtual ~SphericalOffset();

    static SphericalOffsetSPtr create(double offset);

    double getOffset() const;
    void setOffset(double offset);
    bool isInfJump() const;
    void setInfJump(bool inf_jump);

    int compareTo(SphericalOffsetSPtr other);

protected:
    SphericalOffset(double offset);
    double offset_;
    bool inf_jump_;
};

} } }

#endif /* DATA_3D_SKEL_SPHERICALOFFSET_H */

