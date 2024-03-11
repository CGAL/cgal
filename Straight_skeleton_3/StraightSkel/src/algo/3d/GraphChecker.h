/**
 * @file   algo/3d/GraphChecker.h
 * @author Gernot Walzl
 * @date   2015-11-04
 */

#ifndef ALGO_3D_GRAPHCHECKER_H
#define ALGO_3D_GRAPHCHECKER_H

#include "algo/3d/ptrs.h"
#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include <set>

namespace algo { namespace _3d {

using namespace data::_3d;
using namespace data::_3d::skel;

class GraphChecker {
public:
    virtual ~GraphChecker();

    static GraphCheckerSPtr create();

    NodeSPtr getNode(AbstractEventSPtr event);
    AbstractEventSPtr findEvent(StraightSkeletonSPtr skel, NodeSPtr node);
    unsigned int countVisitedChilds(const std::set<NodeSPtr>& visited, NodeSPtr node);

    bool check(StraightSkeletonSPtr skel);

protected:
    GraphChecker();
};

} }

#endif /* ALGO_3D_GRAPHCHECKER_H */
