/**
 * @file   db/3d/NodeDAO.h
 * @author Gernot Walzl
 * @date   2013-06-03
 */

#ifndef DB_3D_NODEDAO_H
#define DB_3D_NODEDAO_H

#include "data/3d/ptrs.h"
#include "data/3d/skel/ptrs.h"
#include "data/3d/skel/Node.h"
#include "db/ptrs.h"
#include "db/3d/ptrs.h"
#include "db/3d/DAOFactory.h"
#include <string>

namespace db { namespace _3d {

using data::_3d::skel::Node;
using data::_3d::skel::NodeSPtr;

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

#endif /* DB_3D_NODEDAO_H */
