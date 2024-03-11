#ifndef DB_3D_EDGEDAO_H
#define DB_3D_EDGEDAO_H

#include "data/3d/ptrs.h"
#include "data/3d/Edge.h"
#include "db/ptrs.h"
#include <string>

namespace db { namespace _3d {

using data::_3d::Edge;
using data::_3d::EdgeSPtr;

class EdgeDAO {
friend class DAOFactory;
public:
    virtual ~EdgeDAO();
    std::string getTableSchema() const;
    int insert(EdgeSPtr edge);
    bool del(EdgeSPtr edge);
    EdgeSPtr find(int polyhedronid, int eid);
    bool update(EdgeSPtr edge);
protected:
    EdgeDAO();
    int nextEID(int polyhedronid);
};

} }

#endif /* DB_3D_EDGEDAO_H */

