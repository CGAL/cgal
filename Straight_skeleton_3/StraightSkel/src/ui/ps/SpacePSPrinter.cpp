/**
 * @file   ui/ps/SpacePSPrinter.cpp
 * @author Gernot Walzl
 * @date   2012-11-14
 */

#include "ui/ps/SpacePSPrinter.h"

#include "debug.h"
#include "typedefs_thread.h"
#include "data/3d/KernelFactory.h"
#include "data/3d/Polyhedron.h"
#include "data/3d/Edge.h"
#include "data/3d/Vertex.h"
#include "data/3d/Facet.h"
#include "data/3d/SphericalPolygon.h"
#include "data/3d/CircularEdge.h"
#include "data/3d/CircularVertex.h"
#include "data/3d/skel/StraightSkeleton.h"
#include "data/3d/skel/Arc.h"
#include "data/3d/skel/Node.h"
#include "data/3d/skel/SphericalSkeleton.h"
#include "data/3d/skel/CircularArc.h"
#include "data/3d/skel/CircularNode.h"

namespace ui { namespace ps {

SpacePSPrinter::SpacePSPrinter() {
    // identity matrix
    for (unsigned int i = 0; i < 16; i++) {
        modelview_[i] = 0.0f;
        projection_[i] = 0.0f;
    }
    for (unsigned int i = 0; i < 4; i++) {
        modelview_[i+4*i] = 1.0f;
        projection_[i+4*i] = 1.0f;
    }

    viewport_[0] = 0;
    viewport_[1] = 0;
    viewport_[2] = 640;
    viewport_[3] = 480;
    initBoundingBox();
}

SpacePSPrinter::~SpacePSPrinter() {
    // intentionally does nothing
}

SpacePSPrinterSPtr SpacePSPrinter::create() {
    return SpacePSPrinterSPtr(new SpacePSPrinter());
}

void SpacePSPrinter::initBoundingBox() {
    bounding_box_[0] = viewport_[0];
    bounding_box_[1] = viewport_[1];
    bounding_box_[2] = viewport_[0] + viewport_[2];
    bounding_box_[3] = viewport_[1] + viewport_[3];
}

void SpacePSPrinter::setModelviewMatrix(float modelview[16]) {
    for (unsigned int i = 0; i < 16; i++) {
        modelview_[i] = modelview[i];
    }
}

void SpacePSPrinter::setProjectionMatrix(float projection[16]) {
    for (unsigned int i = 0; i < 16; i++) {
        projection_[i] = projection[i];
    }
}

void SpacePSPrinter::setViewport(int viewport[4]) {
    for (unsigned int i = 0; i < 4; i++) {
        viewport_[i] = viewport[i];
    }
    initBoundingBox();
}

void SpacePSPrinter::setCamEye(const vec3f cam_eye) {
    for (unsigned int i = 0; i < 3; i++) {
        cam_eye_[i] = cam_eye[i];
    }
}

void SpacePSPrinter::setCamCenter(const vec3f cam_center) {
    for (unsigned int i = 0; i < 3; i++) {
        cam_center_[i] = cam_center[i];
    }
}

void SpacePSPrinter::printCommentCamera(std::ostream& out) {
    out << "% cam_eye = <"
        << cam_eye_[0] << ", "
        << cam_eye_[1] << ", "
        << cam_eye_[2] << ">" << std::endl;
    out << "% cam_center = <"
        << cam_center_[0] << ", "
        << cam_center_[1] << ", "
        << cam_center_[2] << ">" << std::endl;
    out << std::endl;
}

void SpacePSPrinter::to2d(const vec3f p3, vec2f& p2) {
    vec4f p_in;
    for (unsigned int i = 0; i < 3; i++) {
        p_in[i] = p3[i];
    }
    p_in[3] = 1.0f;

    vec4f p_mv;
    for (unsigned int r = 0; r < 4; r++) {
        p_mv[r] = 0.0f;
        for (unsigned int c = 0; c < 4; c++) {
            p_mv[r] += modelview_[r + 4*c] * p_in[c];
        }
    }

    vec4f p_proj;
    for (unsigned int r = 0; r < 4; r++) {
        p_proj[r] = 0.0f;
        for (unsigned int c = 0; c < 4; c++) {
            p_proj[r] += projection_[r + 4*c] * p_mv[c];
        }
    }

    float w = p_proj[3];
    vec2f p_nd;  // normalized device coordinates
    for (unsigned int i = 0; i < 2; i++) {
        p_nd[i] = p_proj[i]/w;
    }

    p2[0] = ((p_nd[0] + 1.0f) * ((float)viewport_[2]/2.0f)) + (float)viewport_[0];
    p2[1] = ((p_nd[1] + 1.0f) * ((float)viewport_[3]/2.0f)) + (float)viewport_[1];
}

void SpacePSPrinter::printLine3(const vec3f src, const vec3f dst, std::ostream& out) {
    vec2f src2;
    vec2f dst2;
    to2d(src, src2);
    to2d(dst, dst2);
    printLine(src2, dst2, out);
}

void SpacePSPrinter::printPolyhedron(PolyhedronSPtr polyhedron, std::ostream& out) {
    if (!polyhedron) {
        DEBUG_PRINT("Warning: polyhedron is null.");
        return;
    }
    ReadLock l(polyhedron->mutex());
    std::list<data::_3d::EdgeSPtr>::iterator it_e = polyhedron->edges().begin();
    while (it_e != polyhedron->edges().end()) {
        data::_3d::EdgeSPtr edge = *it_e++;
        data::_3d::VertexSPtr vertex_src = edge->getVertexSrc();
        data::_3d::VertexSPtr vertex_dst = edge->getVertexDst();
        vec3f src = {(float)vertex_src->getX(),
                     (float)vertex_src->getY(),
                     (float)vertex_src->getZ()};
        vec3f dst = {(float)vertex_dst->getX(),
                     (float)vertex_dst->getY(),
                     (float)vertex_dst->getZ()};
        printLine3(src, dst, out);
    }
    out << std::endl;
}


std::list<data::_3d::FacetSPtr> SpacePSPrinter::getFacetsToShade(PolyhedronSPtr polyhedron) {
    std::list<data::_3d::FacetSPtr> result;
    std::list<data::_3d::FacetSPtr> temp;
    std::list<data::_3d::FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        data::_3d::FacetSPtr facet = *it_f++;
        data::_3d::Plane3SPtr plane = facet->plane();
        data::_3d::Vector3SPtr normal = data::_3d::KernelFactory::createVector3(plane);
        vec3f norm = {(float)(*normal)[0],
                (float)(*normal)[1],
                (float)(*normal)[2]};
        float d = 0.0f;
#ifdef USE_CGAL
        d = (float)plane->d();
#else
        d = (float)plane->getD();
#endif
        if (scalar(cam_eye_, norm) + d >= 0) {
            temp.push_back(facet);
        }
    }

    // TODO: this does not work fine
    // a hint is better than nothing
    while (temp.size() > 0) {
        float max_dist = 0.0f;
        std::list<data::_3d::FacetSPtr>::iterator it_max_f = temp.begin();
        it_f = temp.begin();
        while (it_f != temp.end()) {
            data::_3d::FacetSPtr facet = *it_f;
            vec3f center = {0.0f, 0.0f, 0.0f};
            std::list<data::_3d::VertexSPtr>::iterator it_v = facet->vertices().begin();
            while (it_v != facet->vertices().end()) {
                data::_3d::VertexSPtr vertex = *it_v++;
                center[0] += (float)vertex->getX();
                center[1] += (float)vertex->getY();
                center[2] += (float)vertex->getZ();
            }
            scale(center, 1.0f/(float)facet->vertices().size());
            vec3f diff;
            for (unsigned int i = 0; i < 3; i++) {
                diff[i] = center[i] - cam_eye_[i];
            }
            float dist = scalar(diff, diff);
            if (dist > max_dist || max_dist == 0.0f) {
                it_max_f = it_f;
                max_dist = dist;
            }
            it_f++;
        }
        result.push_back(*it_max_f);
        temp.erase(it_max_f);
    }
    return result;
}


void SpacePSPrinter::printPolyhedronShade(PolyhedronSPtr polyhedron,
            float min_gray, float max_gray, bool print_edges, std::ostream& out) {
    if (!polyhedron) {
        DEBUG_PRINT("Warning: polyhedron is null.");
        return;
    }
    vec3f dir_cam;
    for (unsigned int i = 0; i < 3; i++) {
        dir_cam[i] = cam_eye_[i] - cam_center_[i];
    }
    ReadLock l(polyhedron->mutex());
    std::list<data::_3d::FacetSPtr> facets = getFacetsToShade(polyhedron);
    std::list<data::_3d::FacetSPtr>::iterator it_f = facets.begin();
    while (it_f != facets.end()) {
        data::_3d::FacetSPtr facet = *it_f++;
        data::_3d::Plane3SPtr plane = facet->plane();
        data::_3d::Vector3SPtr normal = data::_3d::KernelFactory::createVector3(plane);
        vec3f norm = {(float)(*normal)[0],
                (float)(*normal)[1],
                (float)(*normal)[2]};
        float angle_cam_norm = angle(dir_cam, norm);
        float gray = (cosf(angle_cam_norm) * (max_gray-min_gray)) + min_gray;
        setGray(gray, out);
        std::list<data::_3d::VertexSPtr> vertices;
        data::_3d::EdgeSPtr edge = facet->edges().front();
        data::_3d::EdgeSPtr edge_first;
        while (edge != edge_first) {
            if (!edge_first) {
                edge_first = edge;
            }
            vertices.push_back(edge->src(facet));
            edge = edge->next(facet);
        }
        vec2f points[vertices.size()];
        unsigned int i = 0;
        std::list<data::_3d::VertexSPtr>::iterator it_v = vertices.begin();
        while (it_v != vertices.end()) {
            data::_3d::VertexSPtr vertex = *it_v++;
            vec3f point3 = {(float)vertex->getX(),
                     (float)vertex->getY(),
                     (float)vertex->getZ()};
            to2d(point3, points[i]);
            i++;
        }
        printPath(vertices.size(), points, true, true, out);

        if (print_edges) {
            setGray(0.0f, out);
            std::list<data::_3d::EdgeSPtr>::iterator it_e = facet->edges().begin();
            while (it_e != facet->edges().end()) {
                data::_3d::EdgeSPtr edge = *it_e++;
                data::_3d::VertexSPtr vertex_src = edge->getVertexSrc();
                data::_3d::VertexSPtr vertex_dst = edge->getVertexDst();
                vec3f src = {(float)vertex_src->getX(),
                             (float)vertex_src->getY(),
                             (float)vertex_src->getZ()};
                vec3f dst = {(float)vertex_dst->getX(),
                             (float)vertex_dst->getY(),
                             (float)vertex_dst->getZ()};
                printLine3(src, dst, out);
            }
        }
    }
    out << std::endl;
}

void SpacePSPrinter::printSkel(data::_3d::skel::StraightSkeletonSPtr skel, std::ostream& out) {
    ReadLock l(skel->mutex());
    std::list<data::_3d::skel::ArcSPtr>::iterator it_a = skel->arcs().begin();
    while (it_a != skel->arcs().end()) {
        data::_3d::skel::ArcSPtr arc = *it_a++;
        if (!arc->hasNodeDst()) {
            continue;
        }
        data::_3d::skel::NodeSPtr node_src = arc->getNodeSrc();
        data::_3d::skel::NodeSPtr node_dst = arc->getNodeDst();
        vec3f src = {(float)node_src->getX(),
                     (float)node_src->getY(),
                     (float)node_src->getZ()};
        vec3f dst = {(float)node_dst->getX(),
                     (float)node_dst->getY(),
                     (float)node_dst->getZ()};
        printLine3(src, dst, out);
    }
    out << std::endl;
}

float SpacePSPrinter::to2dRadius(const vec3f center, float radius) {
    vec3f dir_cam;
    for (unsigned int i = 0; i < 3; i++) {
        dir_cam[i] = center[i] - cam_eye_[i];
    }
    float dist = length(dir_cam);
    float alpha = asinf(radius/dist);
    float radius_proj = radius/cosf(alpha);
    vec3f up = {0.0f, 0.0f, 1.0f};
    vec3f dir_p;
    cross(up, dir_cam, dir_p);
    normalize(dir_p);
    scale(dir_p, radius_proj);
    vec3f point3;
    for (unsigned int i = 0; i < 3; i++) {
        point3[i] = center[i] + dir_p[i];
    }
    vec2f point2;
    to2d(point3, point2);
    vec2f center2;
    to2d(center, center2);
    float radius2 = 0.0f;
    for (unsigned int i = 0; i < 2; i++) {
        float diff = point2[i] - center2[i];
        radius2 += diff * diff;
    }
    radius2 = sqrtf(radius2);
    return radius2;
}

void SpacePSPrinter::printSphere(const vec3f center, float radius, std::ostream& out) {
    vec2f center2;
    to2d(center, center2);
    float radius2 = to2dRadius(center, radius);
    printCircle(center2, radius2, out);
    //printCircle(center2, 0.5f, out);
    out << std::endl;
}

void SpacePSPrinter::printCircularEdge(const vec3f center, const vec3f axis,
        const vec3f src, const vec3f dst, std::ostream& out) {
    vec3f dir_src;
    vec3f dir_dst;
    for (unsigned int i = 0; i < 3; i++) {
        dir_src[i] = src[i] - center[i];
        dir_dst[i] = dst[i] - center[i];
    }
    float radius = 0.0f;
    for (unsigned int i = 0; i < 3; i++) {
        radius += dir_src[i] * dir_src[i];
    }
    radius = sqrtf(radius);
    normalize(dir_src);
    normalize(dir_dst);

    unsigned int num_lines = 32;  // has to be a number with base 2
    vec3f dir[num_lines + 1];
    copy(dir_src, dir[0]);
    copy(dir_dst, dir[num_lines]);
    unsigned int inc = num_lines/2;
    while (inc > 0) {
        for (unsigned int i = inc; i < num_lines; i += 2*inc) {
            for (unsigned int j = 0; j < 3; j++) {
                dir[i][j] = dir[i-inc][j] + dir[i+inc][j];
            }
            if (i == num_lines/2) {
                if (length(dir[i]) < 0.0005) {
                    // angle between src and dst = M_PI
                    cross(axis, dir[0], dir[i]);
                } else {
                    // fix orientation when angle between src and dst > M_PI
                    vec3f dir_normal;
                    cross(dir[0], dir[i], dir_normal);
                    if (angle(dir_normal, axis) > M_PI/2.0f) {
                        scale(dir[i], -1.0f);
                    }
                }
            }
            normalize(dir[i]);
        }
        inc /= 2;
    }
    for (unsigned int i = 0; i <= num_lines; i++) {
        scale(dir[i], radius);
    }

    vec2f points2[num_lines + 1];
    for (unsigned int i = 0; i <= num_lines; i++) {
        vec3f point;
        for (unsigned int j = 0; j < 3; j++) {
            point[j] = center[j] + dir[i][j];
        }
        to2d(point, points2[i]);
    }
    printPath(num_lines+1, points2, false, false, out);
}

void SpacePSPrinter::printSphericalPolygon(
        SphericalPolygonSPtr sphericalpolygon, std::ostream& out) {
    data::_3d::Sphere3SPtr sphere = sphericalpolygon->getSphere();
    data::_3d::Point3SPtr p_center =
            data::_3d::KernelFactory::createPoint3(sphere);
    vec3f center = {(float)(*p_center)[0],
                    (float)(*p_center)[1],
                    (float)(*p_center)[2]};
    std::list<data::_3d::CircularEdgeSPtr>::iterator it_e = sphericalpolygon->edges().begin();
    while (it_e != sphericalpolygon->edges().end()) {
        data::_3d::CircularEdgeSPtr edge = *it_e++;
        data::_3d::CircularVertexSPtr vertex_src = edge->getVertexSrc();
        data::_3d::CircularVertexSPtr vertex_dst = edge->getVertexDst();
        if (vertex_src->getPoint() == vertex_dst->getPoint()) {
            continue;
        }
        if (!vertex_src->isPointValid() || !vertex_dst->isPointValid()) {
            continue;
        }
        vec3f src = {(float)vertex_src->getX(),
                     (float)vertex_src->getY(),
                     (float)vertex_src->getZ()};
        vec3f dst = {(float)vertex_dst->getX(),
                     (float)vertex_dst->getY(),
                     (float)vertex_dst->getZ()};
        data::_3d::Vector3SPtr normal = data::_3d::KernelFactory::createVector3(
                edge->getSupportingPlane());
        vec3f axis = {(float)(*normal)[0],
                      (float)(*normal)[1],
                      (float)(*normal)[2]};
        printCircularEdge(center, axis, src, dst, out);
    }
    out << std::endl;
}

void SpacePSPrinter::printSphericalIntersections(
        SphericalPolygonSPtr sphericalpolygon, std::ostream& out) {
    data::_3d::Sphere3SPtr sphere = sphericalpolygon->getSphere();
    data::_3d::Point3SPtr p_center =
            data::_3d::KernelFactory::createPoint3(sphere);
    vec3f center = {(float)(*p_center)[0],
                    (float)(*p_center)[1],
                    (float)(*p_center)[2]};
    vec2f center2;
    to2d(center, center2);
    std::list<data::_3d::CircularVertexSPtr>::iterator it_v = sphericalpolygon->vertices().begin();
    while (it_v != sphericalpolygon->vertices().end()) {
        data::_3d::CircularVertexSPtr vertex = *it_v++;
        vec3f point = {(float)vertex->getX(),
                       (float)vertex->getY(),
                       (float)vertex->getZ()};
        vec2f point2;
        to2d(point, point2);
        vec2f dir2;
        float length2 = 0.0f;
        for (unsigned int i = 0; i < 2; i++) {
            dir2[i] = point2[i] - center2[i];
            length2 += dir2[i] * dir2[i];
        }
        length2 = sqrtf(length2);
        for (unsigned int i = 0; i < 2; i++) {
            dir2[i] /= length2;
        }
        vec2f dst2;
        for (unsigned int i = 0; i < 2; i++) {
            // draw intersecting line 50% longer
            dst2[i] = center2[i] + (dir2[i] * length2 * 1.5f);
        }
        printLine(center2, dst2, out);
    }
    out << std::endl;
}

void SpacePSPrinter::printSphericalSkel(
        data::_3d::skel::SphericalSkeletonSPtr sphericalskel, std::ostream& out) {
    data::_3d::Sphere3SPtr sphere = sphericalskel->getSphere();
    data::_3d::Point3SPtr p_center =
            data::_3d::KernelFactory::createPoint3(sphere);
    vec3f center = {(float)(*p_center)[0],
                    (float)(*p_center)[1],
                    (float)(*p_center)[2]};
    std::list<data::_3d::skel::CircularArcSPtr>::iterator it_a = sphericalskel->arcs().begin();
    while (it_a != sphericalskel->arcs().end()) {
        data::_3d::skel::CircularArcSPtr arc = *it_a++;
        if (!arc->hasNodeDst()) {
            continue;
        }
        data::_3d::skel::CircularNodeSPtr node_src = arc->getNodeSrc();
        data::_3d::skel::CircularNodeSPtr node_dst = arc->getNodeDst();
        vec3f src = {(float)node_src->getX(),
                     (float)node_src->getY(),
                     (float)node_src->getZ()};
        vec3f dst = {(float)node_dst->getX(),
                     (float)node_dst->getY(),
                     (float)node_dst->getZ()};
        data::_3d::Vector3SPtr normal = data::_3d::KernelFactory::createVector3(
                arc->getSupportingPlane());
        vec3f axis = {(float)(*normal)[0],
                      (float)(*normal)[1],
                      (float)(*normal)[2]};
        printCircularEdge(center, axis, src, dst, out);
    }
    out << std::endl;
}

} }
