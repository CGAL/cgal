/**
 * @file   db/3d/OBJFile.h
 * @author Gernot Walzl
 * @date   2012-05-15
 */

#ifndef DB_3D_OBJFILE_H
#define DB_3D_OBJFILE_H

#include "db/3d/AbstractFile.h"
#include <string>

namespace db { namespace _3d {

using namespace data::_3d;

/**
 * Wavefront obj file format
 */
class OBJFile : public AbstractFile {
public:
    virtual ~OBJFile();

    static PolyhedronSPtr load(const std::string& filename);

    /**
     * It is impossible for an obj file to store holes inside a facet.
     */
    static bool save(const std::string& filename, PolyhedronSPtr polyhedron);

protected:
    OBJFile();
};

} }

#endif /* DB_3D_OBJFILE_H */
