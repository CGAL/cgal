#ifndef CGAL_SNC_ITERATION_H
#define CGAL_SNC_ITERATION_H

#undef CGAL_forall_iterators
#define CGAL_forall_iterators(x,L)\
for(x = (L).begin(); x != (L).end(); ++x) 

#define CGAL_forall_vertices(x,SNC)\
for(x = (SNC).vertices_begin(); x != (SNC).vertices_end(); ++x) 
#define CGAL_forall_halfedges(x,SNC)\
for(x = (SNC).halfedges_begin(); x != (SNC).halfedges_end(); ++x)
#define CGAL_forall_edges(x,SNC)\
for(x = (SNC).halfedges_begin(); x != (SNC).halfedges_end(); ++x) \
if ( x->is_twin() ) continue; else
#define CGAL_forall_halffacets(x,SNC)\
for(x = (SNC).halffacets_begin(); x != (SNC).halffacets_end(); ++x) 
#define CGAL_forall_facets(x,SNC)\
for(x = (SNC).halffacets_begin(); x != (SNC).halffacets_end(); ++x) \
if ( x->is_twin() ) continue; else
#define CGAL_forall_volumes(x,SNC)\
for(x = (SNC).volumes_begin(); x != (SNC).volumes_end(); ++x) 
#define CGAL_forall_svertices(x,SNC)\
for(x = (SNC).svertices_begin(); x != (SNC).svertices_end(); ++x) 
#define CGAL_forall_shalfloops(x,SNC)\
for(x = (SNC).shalfloops_begin(); x != (SNC).shalfloops_end(); ++x) 
#define CGAL_forall_shalfedges(x,SNC)\
for(x = (SNC).shalfedges_begin(); x != (SNC).shalfedges_end(); ++x) 
#define CGAL_forall_sedges(x,SNC)\
for(x = (SNC).shalfedges_begin(); x != (SNC).shalfedges_end(); ++(++x))
#define CGAL_forall_sfaces(x,SNC)\
for(x = (SNC).sfaces_begin(); x != (SNC).sfaces_end(); ++x)

#define CGAL_forall_facet_cycles_of(x,F)\
for(x = (F)->facet_cycles_begin(); x != (F)->facet_cycles_end(); ++x) 
#define CGAL_forall_shells_of(x,C)\
for(x = (C)->shells_begin(); x != (C)->shells_end(); ++x) 

#define CGAL_forall_svertices_of(x,V)\
for(x = (V)->svertices_begin(); x != (V)->svertices_end(); ++x)
#define CGAL_forall_sedges_of(x,V)\
for(x = (V)->shalfedges_begin(); x != (V)->shalfedges_end(); ++(++x))
#define CGAL_forall_shalfedges_of(x,V)\
for(x = (V)->shalfedges_begin(); x != (V)->shalfedges_end(); ++x)
#define CGAL_forall_sfaces_of(x,V)\
for(x = (V)->sfaces_begin(); x != (V)->sfaces_end(); ++x)

#define CGAL_forall_sface_cycles_of(x,F)\
for(x = (F)->sface_cycles_begin(); x != (F)->sface_cycles_end(); ++x)

#endif //CGAL_SNC_ITERATION_H

