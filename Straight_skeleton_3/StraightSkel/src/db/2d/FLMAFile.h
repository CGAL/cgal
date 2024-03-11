/**
 * @file   db/2d/FLMAFile.h
 * @author Gernot Walzl
 * @date   2013-12-16
 */

#ifndef DB_2D_FLMAFILE_H
#define DB_2D_FLMAFILE_H

#include "db/2d/AbstractFile.h"
#include <string>

namespace db { namespace _2d {

using namespace data::_2d;

/**
 * AVL flma file format
 */
class FLMAFile : public AbstractFile {
public:
    virtual ~FLMAFile();
    static PolygonSPtr load(const std::string& filename);
    static bool save(const std::string& filename, PolygonSPtr polygon);
protected:
    FLMAFile();
};

} }

#endif /* DB_2D_FLMAFILE_H */
