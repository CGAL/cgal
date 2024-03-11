/**
 * @file   algo/3d/AngleVertexSplitter.cpp
 * @author Gernot Walzl
 * @date   2012-10-01
 */

#include "algo/3d/AngleVertexSplitter.h"

#include "debug.h"
#include "algo/3d/KernelWrapper.h"
#include "data/3d/Vertex.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/skel/SkelVertexData.h"

namespace algo { namespace _3d {

AngleVertexSplitter::AngleVertexSplitter() {
    type_ = AbstractVertexSplitter::ANGLE_VERTEX_SPLITTER;
}

AngleVertexSplitter::~AngleVertexSplitter() {
    // intentionally does nothing
}

AngleVertexSplitterSPtr AngleVertexSplitter::create() {
    AngleVertexSplitterSPtr result =
            AngleVertexSplitterSPtr(new AngleVertexSplitter());
    return result;
}

PolyhedronSPtr AngleVertexSplitter::splitVertex(VertexSPtr vertex) {
    DEBUG_PRINT("WARNING: AngleVertexSplitter::splitVertex(...) does not work in every case.");
    PolyhedronSPtr polyhedron = vertex->getPolyhedron();
    if (vertex->degree() <= 3) {
        return polyhedron;
    }
    std::list<VertexSPtr> queue;
    queue.push_back(vertex);
    while (queue.size() > 0) {
        VertexSPtr vertex1 = queue.front();
        queue.pop_front();

        // 1. search for 2 facets that give a new edge
        double angle_min = 2*M_PI;
        FacetSPtr facet_min1;
        FacetSPtr facet_min2;
        std::list<FacetWPtr>::iterator it_f1 = vertex1->facets().begin();
        while (it_f1 != vertex1->facets().end()) {
            FacetSPtr facet1 = FacetSPtr(*it_f1++);
            std::list<FacetWPtr>::iterator it_f2 = vertex1->facets().begin();
            while (it_f2 != vertex1->facets().end()) {
                FacetSPtr facet2 = FacetSPtr(*it_f2++);
                if (facet1 == facet2 ||
                        facet1->prev(vertex1) == facet2 ||
                        facet1->next(vertex1) == facet2) {
                    continue;
                }
                // for reflex and convex vertices
                double angle = M_PI - KernelWrapper::angle(
                        facet1->plane(), facet2->plane());
                if (angle < angle_min) {
                    facet_min1 = facet1;
                    facet_min2 = facet2;
                    angle_min = angle;
                }
            }
        }

        // 2. split vertex
        VertexSPtr vertex2 = vertex1->split(facet_min1, facet_min2);
        if (vertex1->hasData()) {
            SkelVertexDataSPtr data1 = std::dynamic_pointer_cast<SkelVertexData>(
                    vertex1->getData());
            SkelVertexDataSPtr data2 = SkelVertexData::create(vertex2);
            data2->setNode(data1->getNode());
        }

        // 3. add vertices with degree > 3 to queue
        if (vertex1->degree() > 3) {
            queue.push_back(vertex1);
        }
        if (vertex2->degree() > 3) {
            queue.push_back(vertex2);
        }
        assert(polyhedron->isConsistent());
    }
    return polyhedron;
}

} }
