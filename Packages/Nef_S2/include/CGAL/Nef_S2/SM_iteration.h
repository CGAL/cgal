#ifndef CGAL_SM_ITERATION_H
#define CGAL_SM_ITERATION_H

#define CGAL_forall_iterators(x,S)\
for(x = S.begin(); x != S.end(); ++x)

#define CGAL_forall_vertices(x,SM)\
for(x = (SM).vertices_begin(); x != (SM).vertices_end(); ++x) 

#define CGAL_forall_halfedges(x,SM)\
for(x = (SM).halfedges_begin(); x != (SM).halfedges_end(); ++x) 

#define CGAL_forall_halfloops(x,SM)\
for(x = (SM).halfloops_begin(); x != (SM).halfloops_end(); ++x) 

#define CGAL_forall_faces(x,SM)\
for(x = (SM).faces_begin(); x != (SM).faces_end(); ++x) 

#define CGAL_forall_face_cycles_of(x,F)\
for(x = (F)->face_cycles_begin(); x != (F)->face_cycles_end(); ++x) 

#endif //CGAL_SM_ITERATION_H

