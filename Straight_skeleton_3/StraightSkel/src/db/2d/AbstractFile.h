/**
 * @file   db/2d/AbstractFile.h
 * @author Gernot Walzl
 * @date   2013-12-18
 */

#ifndef DB_2D_ABSTRACTFILE_H
#define DB_2D_ABSTRACTFILE_H

#include "data/2d/ptrs.h"

namespace db { namespace _2d {

using namespace data::_2d;

class AbstractFile {
public:
    virtual ~AbstractFile();
    static bool hasCollinearEdges(VertexSPtr vertex, double epsilon);
    static int mergeCollinearEdges(PolygonSPtr polygon, double epsilon);
protected:
    AbstractFile();
};

} }

#endif /* DB_2D_ABSTRACTFILE_H */

