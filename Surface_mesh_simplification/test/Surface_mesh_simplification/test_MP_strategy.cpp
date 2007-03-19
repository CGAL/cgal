#define STRATEGY_POLICIES "MP_strategy.cpp"

#define VISITOR_CLASS "check_audit_edge_collapse_visitor.cpp"

#define CREATE_VISITOR Visitor lVisitor_(aName.substr(0,aName.find_last_of(".")) + "_MP.audit") ; \
                       Visitor* lVisitor = &lVisitor_ ;

#define VISITOR_ARGUMENT .visitor(lVisitor)

#define CGAL_TESTING_SURFACE_MESH_SIMPLIFICATION_USING_EXTENDED_VISITOR

#include "aux_edge_collapse_test_Polyhedron_3.cpp"

int main(int argc, char **argv)
{
  return aux_main(argc, argv);
}
