/**
 * @file   data/3d/FacetData.h
 * @author Gernot Walzl
 * @date   2011-11-26
 */

#ifndef DATA_3D_FACETDATA_H
#define DATA_3D_FACETDATA_H

#include "data/3d/ptrs.h"

namespace data { namespace _3d {

class FacetData {
public:
    virtual ~FacetData();

    static FacetDataSPtr create(FacetSPtr facet);

    FacetSPtr getFacet() const;
    void setFacet(FacetSPtr facet);

    bool isHighlight() const;
    void setHighlight(bool highlight);

protected:
    FacetData();
    FacetWPtr facet_;
    bool highlight_;
};

} }

#endif /* DATA_3D_FACETDATA_H */

