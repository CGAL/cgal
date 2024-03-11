/**
 * @file   db/2d/VertexDAO.h
 * @author Gernot Walzl
 * @date   2012-01-27
 */

#ifndef DB_2D_VERTEXDAO_H
#define DB_2D_VERTEXDAO_H

#include "data/2d/ptrs.h"
#include "data/2d/Vertex.h"
#include "data/2d/Polygon.h"
#include "db/ptrs.h"
#include "db/2d/ptrs.h"
#include "db/2d/DAOFactory.h"
#include <string>

namespace db { namespace _2d {

using data::_2d::Vertex;
using data::_2d::VertexSPtr;

class VertexDAO {
friend class DAOFactory;
public:
    virtual ~VertexDAO();
    std::string getTableSchema() const;
    int insert(VertexSPtr vertex);
    bool del(VertexSPtr vertex);
    VertexSPtr find(int polyid, int vid);
    bool update(VertexSPtr vertex);
protected:
    VertexDAO();
    int nextVID(int polyid);
};

} }

#endif /* DB_2D_VERTEXDAO_H */

