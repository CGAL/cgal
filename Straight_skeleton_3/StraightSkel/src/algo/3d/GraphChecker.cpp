/**
 * @file   algo/3d/GraphChecker.cpp
 * @author Gernot Walzl
 * @date   2015-11-04
 */

#include "algo/3d/GraphChecker.h"

#include "debug.h"
#include "data/3d/Vertex.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/Arc.h"
#include "data/3d/skel/StraightSkeleton.h"
#include "data/3d/skel/AbstractEvent.h"
#include "data/3d/skel/EdgeEvent.h"
#include "data/3d/skel/EdgeMergeEvent.h"
#include "data/3d/skel/TriangleEvent.h"
#include "data/3d/skel/DblEdgeMergeEvent.h"
#include "data/3d/skel/DblTriangleEvent.h"
#include "data/3d/skel/TetrahedronEvent.h"
#include "data/3d/skel/VertexEvent.h"
#include "data/3d/skel/FlipVertexEvent.h"
#include "data/3d/skel/SurfaceEvent.h"
#include "data/3d/skel/PolyhedronSplitEvent.h"
#include "data/3d/skel/SplitMergeEvent.h"
#include "data/3d/skel/EdgeSplitEvent.h"
#include "data/3d/skel/PierceEvent.h"
#include <list>

namespace algo { namespace _3d {

GraphChecker::GraphChecker() {
}

GraphChecker::~GraphChecker() {
}

GraphCheckerSPtr GraphChecker::create() {
    GraphCheckerSPtr result = GraphCheckerSPtr(new GraphChecker());
    return result;
}

NodeSPtr GraphChecker::getNode(AbstractEventSPtr event) {
    NodeSPtr result;
    if (event->getType() == AbstractEvent::EDGE_EVENT) {
        result = std::dynamic_pointer_cast<EdgeEvent>(event)->getNode();
    } else if (event->getType() == AbstractEvent::EDGE_MERGE_EVENT) {
        result = std::dynamic_pointer_cast<EdgeMergeEvent>(event)->getNode();
    } else if (event->getType() == AbstractEvent::TRIANGLE_EVENT) {
        result = std::dynamic_pointer_cast<TriangleEvent>(event)->getNode();
    } else if (event->getType() == AbstractEvent::DBL_EDGE_MERGE_EVENT) {
        result = std::dynamic_pointer_cast<DblEdgeMergeEvent>(event)->getNode();
    } else if (event->getType() == AbstractEvent::DBL_TRIANGLE_EVENT) {
        result = std::dynamic_pointer_cast<DblTriangleEvent>(event)->getNode();
    } else if (event->getType() == AbstractEvent::TETRAHEDRON_EVENT) {
        result = std::dynamic_pointer_cast<TetrahedronEvent>(event)->getNode();
    } else if (event->getType() == AbstractEvent::VERTEX_EVENT) {
        result = std::dynamic_pointer_cast<VertexEvent>(event)->getNode();
    } else if (event->getType() == AbstractEvent::FLIP_VERTEX_EVENT) {
        result = std::dynamic_pointer_cast<FlipVertexEvent>(event)->getNode();
    } else if (event->getType() == AbstractEvent::SURFACE_EVENT) {
        result = std::dynamic_pointer_cast<SurfaceEvent>(event)->getNode();
    } else if (event->getType() == AbstractEvent::POLYHEDRON_SPLIT_EVENT) {
        result = std::dynamic_pointer_cast<PolyhedronSplitEvent>(event)->getNode();
    } else if (event->getType() == AbstractEvent::SPLIT_MERGE_EVENT) {
        result = std::dynamic_pointer_cast<SplitMergeEvent>(event)->getNode();
    } else if (event->getType() == AbstractEvent::EDGE_SPLIT_EVENT) {
        result = std::dynamic_pointer_cast<EdgeSplitEvent>(event)->getNode();
    } else if (event->getType() == AbstractEvent::PIERCE_EVENT) {
        result = std::dynamic_pointer_cast<PierceEvent>(event)->getNode();
    }
    return result;
}

AbstractEventSPtr GraphChecker::findEvent(StraightSkeletonSPtr skel, NodeSPtr node) {
    AbstractEventSPtr result;
    std::list<AbstractEventSPtr>::iterator it_e = skel->events().begin();
    while (it_e != skel->events().end()) {
        AbstractEventSPtr event = *it_e++;
        if (getNode(event) == node) {
            result = event;
            break;
        }
    }
    return result;
}

unsigned int GraphChecker::countVisitedChilds(const std::set<NodeSPtr>& visited, NodeSPtr node) {
    unsigned int result = 0;
    std::list<ArcWPtr>::iterator it_a_wptr = node->arcs().begin();
    while (it_a_wptr != node->arcs().end()) {
        ArcWPtr arc_wptr = *it_a_wptr++;
        if (arc_wptr.expired()) continue;
        ArcSPtr arc(arc_wptr);
        NodeSPtr other;
        if (arc->getNodeSrc() == node) {
            other = arc->getNodeDst();
        } else if (arc->getNodeDst() == node) {
            other = arc->getNodeSrc();
        } else {
            DEBUG_VAR(arc->toString());
        }
        if (visited.find(other) != visited.end()) {
            result++;
        }
    }
    return result;
}


bool GraphChecker::check(StraightSkeletonSPtr skel) {
    bool result = true;

    std::set<NodeSPtr> hidden;
    std::set<NodeSPtr> visited;

    // initialization: visit all nodes on boundary
    PolyhedronSPtr polyhedron = skel->getPolyhedron();
    std::list<NodeSPtr>::iterator lit_n = skel->nodes().begin();
    while (lit_n != skel->nodes().end()) {
        NodeSPtr node = *lit_n++;
        bool is_boundary = false;
        std::list<VertexSPtr>::iterator lit_v = polyhedron->vertices().begin();
        while (lit_v != polyhedron->vertices().end()) {
            VertexSPtr vertex = *lit_v++;
            if (node->getPoint() == vertex->getPoint()) {
                is_boundary = true;
                break;
            }
        }
        if (is_boundary) {
            visited.insert(node);
        } else {
            hidden.insert(node);
        }
    }

    int visit_event_types = 1;    // 1 for vanish events; 2 for contact events
    unsigned int num_layers_vanish = 1;
    unsigned int num_layers_contact = 0;
    double min_offset = 0.0;
    unsigned int num_nodes_layer = 0;
    unsigned int size_hidden_begin = 0;
    while (hidden.size() > 0) {
        size_hidden_begin = hidden.size();
        DEBUG_VAR(hidden.size());
        std::set<NodeSPtr>::iterator it_n = hidden.begin();
        while (it_n != hidden.end()) {
            NodeSPtr node = *it_n++;

            int event_type = findEvent(skel, node)->getType();
            if (visit_event_types == 1 && event_type > 7) {
                // interested in vanish events, but contact event.
                continue;
            } else if (visit_event_types == 2 && event_type <= 7) {
                // interested in contact events, but vanish event.
                continue;
            }

            if (visit_event_types == 2 && node->getOffset() < min_offset) {
                continue;
            }

            unsigned int num_visited_childs = countVisitedChilds(visited, node);
            if (num_visited_childs >= 2) {
                if (visit_event_types == 1) {
                    if (node->getOffset() < min_offset) {
                        min_offset = node->getOffset();
                    }
                }
                num_nodes_layer++;
                visited.insert(node);
                hidden.erase(node);
            }
        }
        if (size_hidden_begin == hidden.size()) {
            if (num_nodes_layer == 0) {
                result = false;
                break;
            }
            if (visit_event_types == 1) {
                DEBUG_VAR(num_nodes_layer);
                DEBUG_PRINT("-- No more vanish events");
                visit_event_types = 2;
                DEBUG_PRINT("-- Visiting contact events");
                num_nodes_layer = 0;
                num_layers_contact++;
                DEBUG_VAR(min_offset);
            } else if (visit_event_types == 2) {
                DEBUG_VAR(num_nodes_layer);
                DEBUG_PRINT("-- No more contact events");
                visit_event_types = 1;
                DEBUG_PRINT("-- Visiting vanish events");
                num_nodes_layer = 0;
                num_layers_vanish++;
            }
        }
    }

    DEBUG_VAR(hidden.size());
    DEBUG_VAR(num_layers_vanish);
    DEBUG_VAR(num_layers_contact);
    DEBUG_VAR(result);

    if (result == 0) {
        std::set<NodeSPtr>::iterator it_n = hidden.begin();
        while (it_n != hidden.end()) {
            NodeSPtr node = *it_n++;
            unsigned int num_visited_childs = countVisitedChilds(visited, node);
            AbstractEventSPtr event = findEvent(skel, node);
            DEBUG_VAR(event->toString());
            DEBUG_VAR(num_visited_childs);
        }
    }

    return result;
}

} }
