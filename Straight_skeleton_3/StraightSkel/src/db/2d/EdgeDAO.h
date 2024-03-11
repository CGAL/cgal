/**
 * @file   db/2d/EdgeDAO.h
 * @author Gernot Walzl
 * @date   2012-01-27
 */

#ifndef DB_2D_EDGEDAO_H
#define DB_2D_EDGEDAO_H

#include "data/2d/ptrs.h"
#include "data/2d/Edge.h"
#include "data/2d/skel/SkelEdgeData.h"
#include "db/ptrs.h"
#include "db/2d/ptrs.h"
#include "db/2d/DAOFactory.h"
#include <string>

namespace db { namespace _2d {

using data::_2d::Edge;
using data::_2d::EdgeSPtr;
using data::_2d::skel::SkelEdgeData;
using data::_2d::skel::SkelEdgeDataSPtr;

class EdgeDAO {
friend class DAOFactory;
public:
    virtual ~EdgeDAO();
    std::string getTableSchema() const;
    std::string getTable2Schema() const;
    int insert(EdgeSPtr edge);
    bool del(EdgeSPtr edge);
    EdgeSPtr find(int polyid, int eid);
    bool update(EdgeSPtr edge);
protected:
    EdgeDAO();
    int nextEID(int polyid);
};

} }

#endif /* DB_2D_EDGEDAO_H */
