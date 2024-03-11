/**
 * @file   ui/ps/PSPrinter.cpp
 * @author Gernot Walzl
 * @date   2012-11-06
 */

#include "ui/ps/PSPrinter.h"

namespace ui { namespace ps {

PSPrinter::PSPrinter() {
    // A4
    bounding_box_[0] = 0;
    bounding_box_[1] = 0;
    bounding_box_[2] = 595;
    bounding_box_[3] = 842;
}

PSPrinter::~PSPrinter() {
    // intentionally does nothing.
}

void PSPrinter::printHead(std::ostream& out) {
    out << "%!PS-Adobe-3.0 EPSF-3.0" << std::endl;
    out << "%%BoundingBox: "
        << bounding_box_[0] << ' ' << bounding_box_[1] << ' '
        << bounding_box_[2] << ' ' << bounding_box_[3] << std::endl;
    out << std::endl;
}

void PSPrinter::printComment(const std::string& comment, std::ostream& out) {
    out << "% " << comment << std::endl;
}

void PSPrinter::setLineWidth(float linewidth, std::ostream& out) {
    out << linewidth << " setlinewidth" << std::endl;
}

void PSPrinter::setGray(float gray, std::ostream& out) {
    out << gray << " setgray" << std::endl;
}

void PSPrinter::printLine(const vec2f src, const vec2f dst, std::ostream& out) {
    out << "newpath" << std::endl;
    out << src[0] << ' ' << src[1] << " moveto" << std::endl;
    out << dst[0] << ' ' << dst[1] << " lineto" << std::endl;
    out << "stroke" << std::endl;
}

void PSPrinter::printPath(unsigned int num_points, const vec2f points[],
            bool closepath, bool fill, std::ostream& out) {
    out << "newpath" << std::endl;
    for (unsigned int i = 0; i < num_points; i++) {
        out << points[i][0] << ' ' << points[i][1];
        if (i == 0) {
            out << " moveto";
        } else {
            out << " lineto";
        }
        out << std::endl;
    }
    if (closepath || fill) {
        out << "closepath" << std::endl;
    }
    if (fill) {
        out << "fill" << std::endl;
    } else {
        out << "stroke" << std::endl;
    }
}

void PSPrinter::printCircle(const vec2f center, float radius, std::ostream& out) {
    out << center[0] << ' ' << center[1] << ' '
        << radius << ' '
        << "0 360 arc" << std::endl;
    out << "stroke" << std::endl;
}

} }
