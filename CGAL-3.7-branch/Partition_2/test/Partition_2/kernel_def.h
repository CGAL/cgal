#ifndef REP_DEF_H
#define REP_DEF_H

#if TESTR==1
   typedef double               NT;
   typedef CGAL::Cartesian<NT>  K;
#endif 

#if TESTR==2
   typedef double                NT;
   typedef CGAL::Homogeneous<NT> K;
#endif


#endif // REP_DEF_H
