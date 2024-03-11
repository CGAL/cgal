/**
 * @file   data/3d/skel/Sheet.h
 * @author Gernot Walzl
 * @date   2012-03-27
 */

#ifndef DATA_3D_SKEL_SHEET_H
#define DATA_3D_SKEL_SHEET_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include <list>
#include <memory>
#include <string>

namespace data { namespace _3d { namespace skel {

class Sheet : public std::enable_shared_from_this<Sheet> {
public:
    virtual ~Sheet();
    static SheetSPtr create();

    FacetSPtr getFacetB() const;
    void setFacetB(FacetSPtr facet_b);
    FacetSPtr getFacetF() const;
    void setFacetF(FacetSPtr facet_f);

    StraightSkeletonSPtr getSkel() const;
    void setSkel(StraightSkeletonSPtr skel);
    std::list<SheetSPtr>::iterator getListIt() const;
    void setListIt(std::list<SheetSPtr>::iterator list_it);

    int getID() const;
    void setID(int id);

    Plane3SPtr getPlane() const;
    void setPlane(Plane3SPtr plane);

    void addNode(NodeSPtr node);
    bool removeNode(NodeSPtr node);

    void addArc(ArcSPtr arc);
    bool removeArc(ArcSPtr arc);

    std::list<ArcSPtr>& arcs();
    std::list<NodeSPtr>& nodes();

    std::string toString() const;

protected:
    Sheet();
    FacetSPtr facet_b_;
    FacetSPtr facet_f_;
    std::list<ArcSPtr> arcs_;
    std::list<NodeSPtr> nodes_;
    StraightSkeletonWPtr skel_;
    std::list<SheetSPtr>::iterator list_it_;
    Plane3SPtr plane_;
    int id_;
};

} } }

#endif /* DATA_3D_SKEL_SHEET_H */
