#ifndef CGAL_SM_ITERATION_H
#define CGAL_SM_ITERATION_H

#undef CGAL_forall_iterators
#define CGAL_forall_iterators(x,S)\
for(x = (S).begin(); x != (S).end(); ++x)

#undef CGAL_forall_svertices
#define CGAL_forall_svertices(x,SM)\
for(x = (SM).svertices_begin(); x != (SM).svertices_end(); ++x) 

#undef CGAL_forall_shalfedges
#define CGAL_forall_shalfedges(x,SM)\
for(x = (SM).shalfedges_begin(); x != (SM).shalfedges_end(); ++x) 

#undef CGAL_forall_sedges
#define CGAL_forall_sedges(x,SM)\
for(x = (SM).shalfedges_begin(); x != (SM).shalfedges_end(); ++(++x))

#undef CGAL_forall_shalfloops
#define CGAL_forall_shalfloops(x,SM)\
for(x = (SM).shalfloops_begin(); x != (SM).shalfloops_end(); ++x) 

#undef CGAL_forall_sfaces
#define CGAL_forall_sfaces(x,SM)\
for(x = (SM).sfaces_begin(); x != (SM).sfaces_end(); ++x) 

#undef CGAL_forall_sface_cycles_of
#define CGAL_forall_sface_cycles_of(x,F)\
for(x = (F)->sface_cycles_begin(); x != (F)->sface_cycles_end(); ++x) 

#endif //CGAL_SM_ITERATION_H

