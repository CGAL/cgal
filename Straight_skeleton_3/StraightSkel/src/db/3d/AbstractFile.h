/**
 * @file   db/3d/AbstractFile.h
 * @author Gernot Walzl
 * @date   2013-12-18
 */

#ifndef DB_3D_ABSTRACTFILE_H
#define DB_3D_ABSTRACTFILE_H

#include "data/3d/ptrs.h"

namespace db { namespace _3d {

using namespace data::_3d;

class AbstractFile {
public:
    virtual ~AbstractFile();

    static bool hasCoplanarFacets(EdgeSPtr edge, double epsilon);
    static int mergeCoplanarFacets(PolyhedronSPtr polyhedron, double epsilon);
    static int removeVerticesDegLt3(PolyhedronSPtr polyhedron);
protected:
    AbstractFile();
};

} }

#endif /* DB_3D_ABSTRACTFILE_H */
