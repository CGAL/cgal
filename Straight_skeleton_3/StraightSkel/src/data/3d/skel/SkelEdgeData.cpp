/**
 * @file   data/3d/skel/SkelEdgeData.cpp
 * @author Gernot Walzl
 * @date   2012-04-05
 */

#include "data/3d/skel/SkelEdgeData.h"

#include "debug.h"
#include "data/3d/Edge.h"

namespace data { namespace _3d { namespace skel {

SkelEdgeData::SkelEdgeData() {
    // intentionally does nothing
}

SkelEdgeData::~SkelEdgeData() {
    // intentionally does nothing
}

SkelEdgeDataSPtr SkelEdgeData::create(EdgeSPtr edge) {
    SkelEdgeDataSPtr result = SkelEdgeDataSPtr(new SkelEdgeData());
    result->setEdge(edge);
    edge->setData(result);
    return result;
}

SheetSPtr SkelEdgeData::getSheet() const {
    // DEBUG_WPTR(sheet_);
    if (this->sheet_.expired())
        return SheetSPtr();
    else
        return SheetSPtr(this->sheet_);
}

void SkelEdgeData::setSheet(SheetSPtr sheet) {
    this->sheet_ = sheet;
}

EdgeSPtr SkelEdgeData::getOffsetEdge() const {
    DEBUG_WPTR(offset_edge_);
    if (this->offset_edge_.expired())
        return EdgeSPtr();
    else
        return EdgeSPtr(this->offset_edge_);
}

void SkelEdgeData::setOffsetEdge(EdgeSPtr offset_edge) {
    this->offset_edge_ = offset_edge;
}

FacetSPtr SkelEdgeData::getFacetOrigin() const {
    DEBUG_WPTR(facet_origin_);
    if (this->facet_origin_.expired())
        return FacetSPtr();
    else
        return FacetSPtr(this->facet_origin_);
}

void SkelEdgeData::setFacetOrigin(FacetSPtr facet_origin) {
    this->facet_origin_ = facet_origin;
}

} } }
