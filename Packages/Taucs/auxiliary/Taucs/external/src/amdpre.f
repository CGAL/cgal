c---------------------------------------------------------
c ftp://ftp.cise.ufl.edu/pub/faculty/davis/AMD/amdpre.f
c---------------------------------------------------------
c
c	AMDPRE:  approximate minimum degree ordering
c	algorithm.  Removes "dense" nodes and then
c	calls AMDBAR.  See the tech report describing this
c	code at:
c
c	ftp://ftp.cise.ufl.edu/pub/faculty/davis/AMD/amdpre.ps
c
c       Written by:  Dr. Tim Davis and Joseph L Carmen.
c	davis@cise.ufl.edu
c
c	The primary purpose of this preprossor program is 
c	to detect dense nodes and partition the matrix into  
c       four quadrants.  Where the top left quadrant holds the 
c       sparse nodes and the bottom right quadrant holds the 
c       dense nodes. The top left is then sent to the AMD program 
c       which returns an ordering.  The AMDpre orders the bottom 
c       right in degree order, and returns the ordering for the 
c       entire matrix.
c
c	May 1, 1997
c
c	NOTE:  This routine calls AMDBAR.  It can easily
c	be modified to call the other AMD routines.
c
c---------------------------------------------------------

	subroutine amdpre
     $		(n, pe, iw, len, iwlen, pfree, nv, next,
     $		last, head, elen, degree, ncmpa, w, iovflo,
     $          mapping)

	integer n, iwlen, pfree, ncmpa, iovflo, iw (iwlen), pe (n),
     $		degree (n), nv (n), next (n), last (n), head (n),
     $		elen (n), w (n), len (n), mapping (n)

c--------------------------------------------------------
c
c
c n:	The matrix order.
c
c
c iwlen: The length of iw (1..iwlen).  On input, the matrix is
c	 stored in iw (1..pfree-1).  However, iw (1..iwlen) should be
c	 slightly larger than what is required to hold the matrix, at
c	 least iwlen .ge. pfree + n is recommended. 
c
c pe:	On input, pe (i) is the index in iw of the start of row i, or
c	zero if row i has no off-diagonal non-zeros.  Must of these
c	values will changed if the iw array is compressed.
c
c	
c pfree:  On input the tail end of the array, iw (pfree..iwlen),
c	  is empty, and the matrix is stored in iw (1..pfree-1).  This 
c	  will change if any rows are removed.
c
c
c
c len:  On input, len (i) holds the number of entries in row i of the
c	matrix, excluding the diagonal.  The contents of len (1..n)
c	are undefined on output.  Some entries will change if rows
c	are removed.

c iw:	On input, iw (1..pfree-1) holds the description of each row i
c	in the matrix.  The matrix must be symmetric, and both upper
c	and lower triangular parts must be present.  The diagonal must
c	not be present.  Row i is held as follows:
c
c		len (i):  the length of the row i data structure
c		iw (pe (i) ... pe (i) + len (i) - 1):
c
c		Note that the rows need not be in any particular order,
c		and there may be empty space between the rows.
c
c last:	On output, last (1..n) holds the permutation (the same as the
c	'PERM' argument in Sparspak).  That is, if i = last (k), then
c	row i is the kth pivot row.  Row last (k) of A is the k-th row
c	in the permuted matrix, PAP^T.
c
c elen:	On output elen (1..n) holds the inverse permutation (the same
c	as the 'INVP' argument in Sparspak).  That is, if k = elen (i),
c	then row i is the kth pivot row.  Row i of A appears as the
c	(elen(i))-th row in the permuted matrix, PAP^T.
c	During execution, elen(i) holds the node in the matrix and   
c	is divided into two parts: 
c
c head:	During execution, head(i) holds the nodes of degree i, where 
c	i > dense  and i <= n.  The only entries in the head are nodes that
c	will be removed from the iw array.  head(i) is the starting point
c	for a linked list to the next(i) pointer array.
c
c next:	During execution, is a linked list where next(i) holds 
c	pointers to next(j) where i != j. If next(i) == 0 then
c	i is the last node in the list which started at head(j). 
c	
c mapping: 	The single most important array in the preprocessor.
c		the mapping array is the inverse of the elen array.
c		This array cannot be changed in the AMD program. The
c		mapping array is used to convert the nodes in the 
c		last(n) array returned from the AMD program to their
c		original value. 
c		(need not be defined by the user on input) 
c
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c
c       Local declarations
c
c       The first row of the integer list is required to be
c       saved through the call to the AMD.  The rest are just
c       control variables
c
c
c---------------------------------------------------------------------

	integer ntemp, flag,
     $          lastnode, dense, current, pnum, deg,
     $          number, node, i, j
	real  Z

c--------------------------------------------------------------
c
c Z:	 The variable Z has two functions:
c        1) When Z is set equal to 0 the preprocessor will be 
c           bypassed.
c        2) Z is also used to adjust dense.  The value given to Z
c	    depends on the matrix and will adjust the value of dense
c           where dense = sqrt(n) * Z.  The default value for Z is 1.0
c           The calling program should be modified to pass in this value.
c
c lastnode: The final value of lastnode is the number of nodes
c           sent to the AMD.
c
c flag:	 initially set equal to 0.  If the preprocessor detects a dense
c        row flag will then be set equal to 1.
c
c ntemp: Before the call to the AMD, ntemp is used to save the original
c	 value of n.
c
c dense: This is the key to the preprocessor.  A good value to dense 
c	 will give good results, however, there is no algorithm that will 
c	 select the optimal dense.  A common dense to choose is where
c	 dense = sqrt(n) * Z. 
c
c
c pnum:       value of previous node in degree linked list, also used
c             as a pointer to entries in the iw array.
c number:     value of node in current position of linked list
c node:       temporary storage of a node
c current:    used to track position in an array
c deg:        temporary storage of a node's degree
c i:          do loop control
c j:          do loop control
c
c----------------------------------------------------------------------- 

c----------------------------------------------------------------------
c
c       User can change the value of Z to adjust dense here
c
c----------------------------------------------------------------------
	
	Z = 1.0

c----------------------------------------------------------------------
c
c       Start AMD preprocessing
c
c----------------------------------------------------------------------

c       ** do not change the value of flag
	flag = 0

	if (Z .gt. 0) then

c----------------------------------------------------------------------
c
c	Compute dense.
c       
c---------------------------------------------------------------------

c             ** set dense equal to sqrt(n)
	      dense = Z * sqrt(real(n)) 

c---------------------------------------------------------------------
c	 initialize head(n) and next(n)
c---------------------------------------------------------------------

	   do 10 i = 1, n
	      head(i) = 0
	      next(i) = 0
 10	   continue
	   

c---------------------------------------------------------------------
c        create the degree hash buckets and linked lists 
c        for the dense nodes
c---------------------------------------------------------------------

	   do 20 i = 1, n
	      deg = len(i)
	      if ( deg  .gt. dense) then
		 
c                ** a dense row was found 
		 flag = 1

c                ** insert node in degree list
		 next(i) = head(deg)
		 head(deg) = i

	      end if

 20	   continue

c---------------------------------------------------------------------
c
c       1) Recalculate the degree length of all nodes adjacent to
c       the dense nodes in the degree list.  (Note:  Many of the 
c       dense nodes in the degree list will no longer be dense after
c       this section.)
c
c       2) Constuct the ordering for the nodes not sent to AMD by 
c       selecting the most dense node in the degree list and 
c       then reduce the lengths of all adjacent nodes. Repeat this
c       until no nodes are left with length higher than dense.
c       The dense nodes are placed in the last(n) array.
c          NOTE:  1) nodes are placed after the final value
c                     of lastnode in the last(n) array
c                 2) the AMD routine will not effect anything after lastnode
c                    in the last(n) array.
c                 3) nodes are saved in degree order and in thier original 
c                    state, i.e., no reverse mapping is needed on these.
c---------------------------------------------------------------------


	   if (flag .eq. 1) then
	      
	      lastnode = n
	      dense = dense + 1
	      current = n
	      
c             ** get node from bucket
 40	      node = head(current)


c             ** main loop control
 60	      if (node .eq. 0) then
		 current = current - 1
		 if (current .lt. dense) then
		    go to 70
		 else
		    go to 40
		 endif
	      endif

c   	      ** remove node from bucket
	      head(current) = next(node)
		 
c             ** get degree of current node
	      deg = len(node)


		 
c             ** skip this node if degree was changed to less than dense
	      if (deg .lt. dense) then
		 go to 40
	      endif
		 
c             ** check if degree was changed
	      if (deg .lt. current) then
		 
c                ** insert back into linked list at the lower degree
		 next(node) = head(deg)
		 head(deg) = node
		 
	      else
		    
c                ** insert into last(n)
		 last(lastnode) = node
		 lastnode = lastnode - 1
		    
c                ** len is flagged for use in the mapping contruction 
		 len(node) = 2 * n
		    
c                ** update degree lengths of adjacent nodes
		 if (node .lt. n) then
		    pnum = pe(node + 1) - 1
		 else
		    pnum = pfree - 1
		 endif
		 do 65 i = pe(node), pnum
		    number = iw(i)
		    len(number) = len(number) - 1
 65		 continue
	      endif


	      go to 40

 70	      continue
	   
c---------------------------------------------------------------------
c ************  begin loop to contruct the mapping array
c                the mapping array will place the low dense nodes
c                at the begining and the high dense rows at the end
c                the mapping array is basically a renumbering of the 
c                nodes.
c       ***  NOTE:
c                 forward mapping == elen(n)
c                 reverse mapping == mapping(n)
c---------------------------------------------------------------------
	      lastnode = n
	      current = 1
	      
	      do 80 i = 1, n
		 deg = len(i)
		 if (deg .lt. dense) then
		    
c                   ** insert node at beginning part of elen array
		    elen(i) = current
		    mapping(current) = i
		    current = current + 1
		 else
		    
c                   ** insert node at end part of elen array
		    elen(i) = lastnode
		    mapping(lastnode) = i
		    lastnode = lastnode - 1
		 endif
 80	      continue
	      
c---------------------------------------------------------------------
c *********  construct the new iw array 
c       include only the nodes that are less than or equal to
c       lastnode in the iw array.  lastnode is currently
c       equal to the highest node value that will go to
c       the amd routine.  elen is used for the forward mapping.
c---------------------------------------------------------------------

	      current = 1
	      node = 1
	      
	      do 90 i = 1 , n-1
		 
c                ** compare forward mapping on node i to lastnode
		 if (elen(i) .le. lastnode) then
		    
c                   **  place node in the new iw array
		    pnum = pe(i)
		    pe(node) = current
		    do 100 j = pnum, pe(i+1)-1 
		       number = elen(iw(j))
		       
c                      ** remove adjacent nodes greater than lastnode
		       if (number .le. lastnode) then
			  iw(current) = number
			  current = current + 1
		       end if 
 100		    continue
		    
c                   ** insert new length of node in len array
		    len(node) = current - pe(node)
		    node = node + 1
		 end if
 90	      continue
	      
c             ** repeat above process for the last node
	      if (elen(n) .le. lastnode) then
		 pnum = pe(n)
		 pe(node) = current
		 do 110 j = pnum, pfree-1
		    number = elen(iw(j))
		    if (number .le. lastnode) then
		       iw(current) = number
		       current = current + 1
		    end if 
 110		 continue
		 len(node) = current - pe(node)
		 node = node + 1
	      end if
	      
	      ntemp = n
	      pfree = current
	      n = lastnode
	   end if
	endif

c---------------------------------------------------------------------
c
c       Call the AMD ordering program
c
c---------------------------------------------------------------------


	call amdbar
     $         (n, pe, iw, len, iwlen, pfree, nv, next,
     $		last, head, elen, degree, ncmpa, w, iovflo)


	
	if (flag .eq. 1) then
	   lastnode = n
	   n = ntemp

c---------------------------------------------------------------------
c        Change nodes in last(1 ... lastnode) to original nodes
c---------------------------------------------------------------------

	   do 120 i = 1, lastnode
	      last(i) = mapping(last(i))
 120	   continue

c---------------------------------------------------------------------
c        Invert last(1 ... n) to elen(1 ... n)
c---------------------------------------------------------------------

	   do 130 i = 1, n
	      number = last(i)
	      elen(number) = i
 130	   continue
	   
	   
	end if


	return
	end
c---------------------------------------------------------------------
c	end of preprocessor subroutine
c---------------------------------------------------------------------





