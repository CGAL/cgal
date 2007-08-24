#include<CGAL/Visibility_complex_2/Point_traits.h>
#include<CGAL/Cartesian.h>
typedef CGAL::Cartesian<int> K;
typedef CGAL::Visibility_complex_2_point_traits<K> Gt;

char * input_disks="\
164 318 \
521 364 \
321 241 \
214 321 \
349 327 \
425 443 \
620 332 \
467 201 \
134 94 \
182 278 \
264 81 \
296 210 \
184 199 \
361 111 \
449 300 \
656 187 \
551 112 \
565 260 \
696 241 \
692 399 \
541 457 \
307 411 \
325 455 \
176 531 \
120 436 \
166 422 \
278 354 \
471 539 \
682 490 \
595 421 \
";

char * input_constraints="\
LR 26 21 \
LR 26 20 \
LL 20 27 \
LL 12 4 \
LR 4 26 \
RR 4 7 \
LR 16 7 \
LL 29 17 \
LL 25 9 \
";

char * path="\
LR 0 25 \
RR 25 27 \
RR 27 20 \
RR 20 1 \
";

int nbit=1740;
int ncbit=481;
