#ifndef REP_DEF_H
#define REP_DEF_H

#if TESTR==1
   typedef double               NT;
   typedef CGAL::Simple_cartesian<NT>  K;
#endif

#if TESTR==2
   typedef double                NT;
   typedef CGAL::Simple_homogeneous<NT> K;
#endif


#endif // REP_DEF_H
