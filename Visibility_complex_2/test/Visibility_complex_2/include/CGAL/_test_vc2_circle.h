#include<CGAL/Visibility_complex_2/Circle_traits.h>
#include<CGAL/Cartesian.h>
#include<CGAL/Gmpz.h>
typedef CGAL::Cartesian<CGAL::Gmpz> K;
typedef CGAL::Visibility_complex_2_circle_traits<K> Gt;

char * input_disks="\
163 443 13000 1 \
422 232 416 1 \
332 313 9092 1 \
534 387 4388 1 \
615 160 7957 1 \
633 296 1000 1 \
402 511 916 1 \
185 64 9748 1 \
412 92 3620 1 \
";

char * input_constraints="\
LL 8 1 \
LR 1 2 \
LL 2 1 \
RR 7 2 \
LL 2 0 \
RR 1 4 \
LL 5 3 \
";

int nbit=83;
int ncbit=43;
