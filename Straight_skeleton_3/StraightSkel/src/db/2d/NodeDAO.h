/**
 * @file   db/2d/NodeDAO.h
 * @author Gernot Walzl
 * @date   2013-05-23
 */

#ifndef DB_2D_NODEDAO_H
#define DB_2D_NODEDAO_H

#include "data/2d/ptrs.h"
#include "data/2d/skel/ptrs.h"
#include "data/2d/skel/Node.h"
#include "db/ptrs.h"
#include "db/2d/ptrs.h"
#include "db/2d/DAOFactory.h"
#include <string>

namespace db { namespace _2d {

using data::_2d::skel::Node;
using data::_2d::skel::NodeSPtr;

class NodeDAO {
friend class DAOFactory;
public:
    virtual ~NodeDAO();
    std::string getTableSchema() const;
    int insert(NodeSPtr node);
    bool del(NodeSPtr node);
    NodeSPtr find(int skelid, int nid);
    bool update(NodeSPtr node);
protected:
    NodeDAO();
    int nextNID(int skelid);
};

} }

#endif /* DB_2D_NODEDAO_H */

