/**
 * @file   ui/ps/CutPatternPrinter.cpp
 * @author Gernot Walzl
 * @date   2012-11-14
 */

#include "ui/ps/CutPatternPrinter.h"

#include "typedefs_thread.h"
#include "data/2d/Polygon.h"
#include "data/3d/Facet.h"
#include "data/3d/Polyhedron.h"
#include "ui/ps/PlanePSPrinter.h"
#include "util/Configuration.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include <list>

namespace ui { namespace ps {

CutPatternPrinter::CutPatternPrinter() {
    // intentionally does nothing
}

CutPatternPrinter::~CutPatternPrinter() {
    // intentionally does nothing
}

void CutPatternPrinter::printCutPattern(PolyhedronSPtr polyhedron, const std::string& filename_prefix) {
    PlanePSPrinterSPtr printer2 = PlanePSPrinter::create();
    double scale = util::Configuration::getInstance()->getDouble(
            "ui_ps_CutPatternPrinter", "scale");
    if (scale > 0.0) {
        printer2->setScale((float)scale);
    }
    ReadLock l(polyhedron->mutex());
    unsigned int i = 0;
    std::list<data::_3d::FacetSPtr>::iterator it_f = polyhedron->facets().begin();
    while (it_f != polyhedron->facets().end()) {
        data::_3d::FacetSPtr facet = *it_f++;
        data::_2d::PolygonSPtr polygon2 = facet->toPolygon();
        printer2->initBoundingBox(polygon2);
        std::stringstream filename_eps;
        filename_eps << filename_prefix << '_'
                << std::setfill('0') << std::setw(3) << i++ << ".eps";
        std::ofstream of;
        of.open(filename_eps.str().c_str());
        if (of.is_open()) {
            printer2->printHead(of);
            printer2->printPolygon(polygon2, of);
            of.close();
        }
    }
}

} }
