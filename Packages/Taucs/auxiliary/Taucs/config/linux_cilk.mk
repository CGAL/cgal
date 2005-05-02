#
# Use Cilk
#
##CC        = $(HOME)/Projects/cilk-devel/current.Bug18/support/cilk-local
##CC        = $(HOME)/Projects/cilk-devel/current/support/cilkclocal
##CFLAGS    = -g -O3 -D_POSIX_C_SOURCE=199506L -x cilk 
##CFLAGS    = -O3 
##LD        = $(CC)
##LDFLAGS   = 

#CILK2C     = $(HOME)/Projects/cilk-devel/current.Bug18/cilk2c/cilk2c 
#CILKINCS   = -I $(HOME)/Projects/cilk-devel/current.Bug18/runtime
#CFLAGS    = -g -O3 

FOUTFLG   = -o ./
COUTFLG   = -o ./

CILKC     = $(HOME)/Projects/cilk-devel/current/support/cilkclocal
CILKFLAGS = -O3 -x cilk
LD        = $(CILKC)
LDFLAGS   = 








