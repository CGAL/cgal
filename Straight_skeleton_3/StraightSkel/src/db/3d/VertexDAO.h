#ifndef DB_3D_VERTEXDAO_H
#define DB_3D_VERTEXDAO_H

#include "data/3d/ptrs.h"
#include "data/3d/Vertex.h"
#include "data/3d/Polyhedron.h"
#include "db/ptrs.h"
#include "db/3d/ptrs.h"
#include "db/3d/DAOFactory.h"
#include <string>

namespace db { namespace _3d {

using data::_3d::Vertex;
using data::_3d::VertexSPtr;

class VertexDAO {
friend class DAOFactory;
public:
    virtual ~VertexDAO();
    std::string getTableSchema() const;
    int insert(VertexSPtr vertex);
    bool del(VertexSPtr vertex);
    VertexSPtr find(int polyhedronid, int vid);
    bool update(VertexSPtr vertex);
protected:
    VertexDAO();
    int nextVID(int polyhedronid);
};

} }

#endif /* DB_3D_VERTEXDAO_H */
