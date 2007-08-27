#include<CGAL/Visibility_complex_2/Polygon_traits.h>
#include<CGAL/Cartesian.h>
typedef CGAL::Cartesian<int> K;
typedef CGAL::Visibility_complex_2_polygon_traits<K> Gt;

const char * input_disks="\
1 10 444 \
1 430 226 \
4 321 482 234 319 347 242 347 342 \
5 96 534 503 467 567 528 487 572 282 592 \
5 299 193 288 149 305 113 355 114 351 152 \
5 423 98 383 56 429 24 485 56 465 84 \
5 46 413 52 325 98 248 182 367 78 474 \
5 232 438 216 410 240 386 262 418 266 442 \
5 674 384 447 377 620 232 696 218 752 262 \
4 451 340 441 320 465 306 475 322 \
4 503 299 493 283 499 267 531 267 \
3 405 299 465 210 467 260 \
3 435 155 692 50 718 113 \
3 190 265 130 112 278 34 \
";

const char * input_constraints="\
LR 2 11 \
RR 9 8 \
LR 10 9 \
LR 5 2 \
RR 2 7 \
RR 6 3 \
";

const char * path="\
LL 0 6 \
LR 6 13 \
RL 13 4 \
LL 4 5 \
LR 5 12 \
RL 12 1 \
";

int nbit=137;
int ncbit=85;
