/**
 * @file   data/3d/skel/SkelFacetData.h
 * @author Gernot Walzl
 * @date   2102-05-07
 */

#ifndef DATA_3D_SKEL_SKELFACETDATA_H
#define DATA_3D_SKEL_SKELFACETDATA_H

#include "data/3d/ptrs.h"
#include "data/3d/FacetData.h"
#include "data/3d/skel/ptrs.h"

namespace data { namespace _3d { namespace skel {

class SkelFacetData : public FacetData {
public:
    virtual ~SkelFacetData();

    static SkelFacetDataSPtr create(FacetSPtr facet);

    FacetSPtr getOffsetFacet() const;
    void setOffsetFacet(FacetSPtr offset_facet);

    FacetSPtr getFacetOrigin() const;
    void setFacetOrigin(FacetSPtr facet_origin);

    double getSpeed() const;
    void setSpeed(double speed);

protected:
    SkelFacetData();
    FacetWPtr offset_facet_;
    FacetWPtr facet_origin_;
    double speed_;
};

} } }

#endif /* DATA_3D_SKEL_SKELFACETDATA_H */
