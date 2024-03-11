/**
 * @file   algo/3d/LineInFacet.cpp
 * @author Gernot Walzl
 * @date   2012-04-26
 */

#include "algo/3d/LineInFacet.h"

#include "data/3d/KernelFactory.h"
#include "data/3d/Facet.h"
#include "data/3d/Edge.h"
#include "data/3d/Vertex.h"
#include "algo/3d/KernelWrapper.h"
#include <list>

namespace algo { namespace _3d {

bool IsLineInFacet(FacetSPtr facet, Line3SPtr line) {
    bool result = false;
    unsigned int crossings = 0;

    Plane3SPtr plane = facet->plane();
    Point3SPtr p_line = KernelFactory::createPoint3(line->point());
    Vector3SPtr v_dir = KernelFactory::createVector3(line);
    Point3SPtr p_line_dir = KernelFactory::createPoint3(line->point() + *v_dir);
    Vector3SPtr v_normal = KernelFactory::createVector3(plane);
    Plane3SPtr plane_ray = KernelFactory::createPlane3(p_line, p_line_dir,
            KernelFactory::createPoint3(line->point() + *v_normal));
    Vector3SPtr v_normal_ray = KernelFactory::createVector3(plane_ray);
    Plane3SPtr plane_orient= KernelFactory::createPlane3(p_line, p_line_dir,
            KernelFactory::createPoint3(line->point() + *v_normal_ray));

    std::list<EdgeSPtr>::iterator it_e = facet->edges().begin();
    while (it_e != facet->edges().end()) {
        EdgeSPtr edge = *it_e++;
        Point3SPtr p_src = edge->getVertexSrc()->getPoint();
        Point3SPtr p_dst = edge->getVertexDst()->getPoint();
        int side_src = KernelWrapper::side(plane_ray, p_src);
        int side_dst = KernelWrapper::side(plane_ray, p_dst);
        if ((side_src > 0 && side_dst < 0) ||
                (side_src < 0 && side_dst > 0)) {
            Point3SPtr p_intersect =
                    KernelWrapper::intersection(plane_ray, edge->line());
            if (p_intersect) {
                if (KernelWrapper::side(plane_orient, p_intersect) > 0) {
                    crossings += 1;
                }
            }
        }
    }

    result = (crossings%2 == 1);
    return result;
}

} }
