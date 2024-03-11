/**
 * @file   ui/ps/CutPatternPrinter.h
 * @author Gernot Walzl
 * @date   2012-11-14
 */

#ifndef UI_PS_CUTPATTERNPRINTER_H
#define UI_PS_CUTPATTERNPRINTER_H

#include "data/2d/ptrs.h"
#include "data/3d/ptrs.h"
#include "ui/ps/ptrs.h"
#include <string>

namespace ui { namespace ps {

using data::_3d::PolyhedronSPtr;

class CutPatternPrinter {
public:
    virtual ~CutPatternPrinter();

    static void printCutPattern(PolyhedronSPtr polyhedron, const std::string& filename_prefix);

protected:
    CutPatternPrinter();
};

} }

#endif /* UI_PS_CUTPATTERNPRINTER_H */
