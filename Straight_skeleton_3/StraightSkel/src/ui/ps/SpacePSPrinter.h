/**
 * @file   ui/ps/SpacePSPrinter.h
 * @author Gernot Walzl
 * @date   2012-11-14
 */

#ifndef UI_PS_SPACEPSPRINTER_H
#define UI_PS_SPACEPSPRINTER_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "ui/vecmath.h"
#include "ui/ps/ptrs.h"
#include "ui/ps/PSPrinter.h"
#include <iostream>
#include <list>

namespace ui { namespace ps {

using data::_3d::PolyhedronSPtr;
using data::_3d::SphericalPolygonSPtr;

class SpacePSPrinter : public PSPrinter {
public:
    virtual ~SpacePSPrinter();

    static SpacePSPrinterSPtr create();

    void initBoundingBox();

    void setModelviewMatrix(float modelview[16]);
    void setProjectionMatrix(float projection[16]);
    void setViewport(int viewport[4]);

    void setCamEye(const vec3f cam_eye);
    void setCamCenter(const vec3f cam_center);

    void printCommentCamera(std::ostream& out);

    /**
     * p_2 = P * V * p_3
     */
    void to2d(const vec3f p3, vec2f& p2);

    void printLine3(const vec3f src, const vec3f dst, std::ostream& out);

    void printPolyhedron(PolyhedronSPtr polyhedron, std::ostream& out);
    void printPolyhedronShade(PolyhedronSPtr polyhedron,
            float min_gray, float max_gray, bool print_edges, std::ostream& out);
    void printSkel(data::_3d::skel::StraightSkeletonSPtr skel, std::ostream& out);

    float to2dRadius(const vec3f center, float radius);

    void printSphere(const vec3f center, float radius, std::ostream& out);
    void printCircularEdge(const vec3f center, const vec3f axis,
            const vec3f src, const vec3f dst, std::ostream& out);

    void printSphericalPolygon(
            SphericalPolygonSPtr sphericalpolygon, std::ostream& out);
    void printSphericalIntersections(
            SphericalPolygonSPtr sphericalpolygon, std::ostream& out);
    void printSphericalSkel(
            data::_3d::skel::SphericalSkeletonSPtr sphericalskel, std::ostream& out);

protected:
    SpacePSPrinter();

    std::list<data::_3d::FacetSPtr> getFacetsToShade(PolyhedronSPtr polyhedron);

    /**
     * index = row + 4*column
     */
    float modelview_[16];
    float projection_[16];
    int viewport_[4];

    vec3f cam_eye_;
    vec3f cam_center_;
};

} }

#endif /* UI_PS_SPACEPSPRINTER_H */
