#include<CGAL/Visibility_complex_2/Ellipse_traits.h>
#include<CGAL/Cartesian.h>
#include<CGAL/Gmpz.h>
typedef CGAL::Cartesian<CGAL::Gmpz> K;
typedef CGAL::Visibility_complex_2_ellipse_traits<K> Gt;

const char * input_disks="\
71565 85244 -117830 -10657967 -4306392 1000000000 \
6330 3390 7282 -4488547 -3634529 1000000000 \
7342 5626 915 -3621505 -3752300 1000000000 \
1182 2900 2257 -1814289 -3245118 1000000000 \
3058 16561 -2745 -3070987 -3022548 1000000000 \
2497 3662 -701 -2643269 -1724262 999999999 \
2911 2166 2363 -3228128 -2160654 1000000000 \
";

const char * input_constraints="\
LR 0 6 \
RL 6 5 \
RR 5 3 \
LL 4 5 \
LR 2 3 \
";

int nbit=58;
int ncbit=33;
