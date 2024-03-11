/**
 * @file   ui/ps/PlanePSPrinter.h
 * @author Gernot Walzl
 * @date   2012-11-14
 */

#ifndef UI_PS_PLANEPSPRINTER_H
#define UI_PS_PLANEPSPRINTER_H

#include "ui/ps/ptrs.h"
#include "ui/ps/PSPrinter.h"
#include "data/2d/ptrs.h"
#include "data/2d/skel/ptrs.h"
#include "data/2d/mesh/ptrs.h"
#include <iostream>

namespace ui { namespace ps {

using data::_2d::PolygonSPtr;

class PlanePSPrinter : public PSPrinter {
public:
    virtual ~PlanePSPrinter();

    static PlanePSPrinterSPtr create();

    float getScale() const;
    void setScale(float scale);

    void initBoundingBox(PolygonSPtr polygon);

    void toPaper(const vec2f in, vec2f& out);

    void printPolygon(PolygonSPtr polygon, std::ostream& out);
    void printSkel(data::_2d::skel::StraightSkeletonSPtr skel, std::ostream& out);
    void printMesh(data::_2d::mesh::MeshSPtr mesh, std::ostream& out);

protected:
    PlanePSPrinter();

    void boundingBoxMin(PolygonSPtr polygon, vec2f& out);
    void boundingBoxMax(PolygonSPtr polygon, vec2f& out);

    vec2f translate_;
    float scale_;
};

} }

#endif /* UI_PS_PLANEPSPRINTER_H */
