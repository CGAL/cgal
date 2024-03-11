/**
 * @file   data/3d/skel/SkelEdgeData.h
 * @author Gernot Walzl
 * @date   2012-04-05
 */

#ifndef DATA_3D_SKEL_SKELEDGEDATA_H
#define DATA_3D_SKEL_SKELEDGEDATA_H

#include "data/3d/ptrs.h"
#include "data/3d/EdgeData.h"
#include "data/3d/skel/ptrs.h"

namespace data { namespace _3d { namespace skel {

class SkelEdgeData : public EdgeData {
public:
    virtual ~SkelEdgeData();

    static SkelEdgeDataSPtr create(EdgeSPtr edge);

    SheetSPtr getSheet() const;
    void setSheet(SheetSPtr sheet);
    EdgeSPtr getOffsetEdge() const;
    void setOffsetEdge(EdgeSPtr offset_edge);

    /**
     * Used by WeightVertexSplitter (intersection with a plane)
     */
    FacetSPtr getFacetOrigin() const;
    void setFacetOrigin(FacetSPtr facet_origin);

protected:
    SkelEdgeData();
    SheetWPtr sheet_;
    EdgeWPtr offset_edge_;
    FacetWPtr facet_origin_;
};

} } }

#endif /* DATA_3D_SKEL_SKELEDGEDATA_H */

