/**
 * @file   ui/ps/PSPrinter.h
 * @author Gernot Walzl
 * @date   2012-11-06
 */

#ifndef UI_PS_PSPRINTER_H
#define UI_PS_PSPRINTER_H

#include "ui/typedefs.h"
#include <iostream>
#include <string>

namespace ui { namespace ps {

class PSPrinter {
public:
    virtual ~PSPrinter();

    void printHead(std::ostream& out);
    void printComment(const std::string& comment, std::ostream& out);
    void setLineWidth(float linewidth, std::ostream& out);
    void setGray(float gray, std::ostream& out);
    void printLine(const vec2f src, const vec2f dst, std::ostream& out);
    void printPath(unsigned int num_points, const vec2f points[],
            bool closepath, bool fill, std::ostream& out);
    void printCircle(const vec2f center, float radius, std::ostream& out);

protected:
    PSPrinter();

    int bounding_box_[4];
};

} }

#endif /* UI_PS_PSPRINTER_H */

