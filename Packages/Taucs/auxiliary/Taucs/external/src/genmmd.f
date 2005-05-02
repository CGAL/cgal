C Sivan: I modified INTEGER*2 -> INTEGER*4

C***************************************************************
C***************************************************************
C****     GENMMD ..... MULTIPLE MINIMUM EXTERNAL DEGREE     ****
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE IMPLEMENTS THE MINIMUM DEGREE
C        ALGORITHM.  IT MAKES USE OF THE IMPLICIT REPRESENTATION
C        OF ELIMINATION GRAPHS BY QUOTIENT GRAPHS, AND THE
C        NOTION OF INDISTINGUISHABLE NODES.  IT ALSO IMPLEMENTS
C        THE MODIFICATIONS BY MULTIPLE ELIMINATION AND MINIMUM
C        EXTERNAL DEGREE.
C        ---------------------------------------------
C        CAUTION - THE ADJACENCY VECTOR ADJNCY WILL BE
C        DESTROYED.
C        ---------------------------------------------
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C        (XADJ,ADJNCY) - THE ADJACENCY STRUCTURE.
C        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
C        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT) INTEGER
C                 (ANY SMALLER ESTIMATE WILL DO) FOR MARKING
C                 NODES.
C
C     OUTPUT PARAMETERS -
C        PERM   - THE MINIMUM DEGREE ORDERING.
C        INVP   - THE INVERSE OF PERM.
C        NOFSUB - AN UPPER BOUND ON THE NUMBER OF NONZERO
C                 SUBSCRIPTS FOR THE COMPRESSED STORAGE SCHEME.
C
C     WORKING PARAMETERS -
C        DHEAD  - VECTOR FOR HEAD OF DEGREE LISTS.
C        INVP   - USED TEMPORARILY FOR DEGREE FORWARD LINK.
C        PERM   - USED TEMPORARILY FOR DEGREE BACKWARD LINK.
C        QSIZE  - VECTOR FOR SIZE OF SUPERNODES.
C        LLIST  - VECTOR FOR TEMPORARY LINKED LISTS.
C        MARKER - A TEMPORARY MARKER VECTOR.
C
C     PROGRAM SUBROUTINES -
C        MMDELM, MMDINT, MMDNUM, MMDUPD.
C
C***************************************************************
C
      SUBROUTINE  GENMMD ( NEQNS, XADJ, ADJNCY, INVP, PERM,
     1                     DELTA, DHEAD, QSIZE, LLIST, MARKER,
     1                     MAXINT, NOFSUB )
C
C***************************************************************
C
C         INTEGER*2  ADJNCY(1), DHEAD(1) , INVP(1)  , LLIST(1) ,
C     1              MARKER(1), PERM(1)  , QSIZE(1)
         INTEGER*4  ADJNCY(1), DHEAD(1) , INVP(1)  , LLIST(1) ,
     1              MARKER(1), PERM(1)  , QSIZE(1)
         INTEGER*4  XADJ(1)
         INTEGER*4  DELTA , EHEAD , I     , MAXINT, MDEG  ,
     1              MDLMT , MDNODE, NEQNS , NEXTMD, NOFSUB,
     1              NUM, TAG
C
C***************************************************************
C
         IF  ( NEQNS .LE. 0 )  RETURN
C
C        ------------------------------------------------
C        INITIALIZATION FOR THE MINIMUM DEGREE ALGORITHM.
C        ------------------------------------------------
         NOFSUB = 0
         CALL  MMDINT ( NEQNS, XADJ, ADJNCY, DHEAD, INVP, PERM,
     1                  QSIZE, LLIST, MARKER )
C
C        ----------------------------------------------
C        NUM COUNTS THE NUMBER OF ORDERED NODES PLUS 1.
C        ----------------------------------------------
         NUM = 1
C
C        -----------------------------
C        ELIMINATE ALL ISOLATED NODES.
C        -----------------------------
         NEXTMD = DHEAD(1)
  100    CONTINUE
             IF  ( NEXTMD .LE. 0 )  GO TO 200
                 MDNODE = NEXTMD
                 NEXTMD = INVP(MDNODE)
                 MARKER(MDNODE) = MAXINT
                 INVP(MDNODE) = - NUM
                 NUM = NUM + 1
                 GO TO 100
C
  200    CONTINUE
C        ----------------------------------------
C        SEARCH FOR NODE OF THE MINIMUM DEGREE.
C        MDEG IS THE CURRENT MINIMUM DEGREE;
C        TAG IS USED TO FACILITATE MARKING NODES.
C        ----------------------------------------
         IF  ( NUM .GT. NEQNS )  GO TO 1000
         TAG = 1
         DHEAD(1) = 0
         MDEG = 2
  300    CONTINUE
             IF  ( DHEAD(MDEG) .GT. 0 )  GO TO 400
                 MDEG = MDEG + 1
                 GO TO 300
  400        CONTINUE
C            -------------------------------------------------
C            USE VALUE OF DELTA TO SET UP MDLMT, WHICH GOVERNS
C            WHEN A DEGREE UPDATE IS TO BE PERFORMED.
C            -------------------------------------------------
             MDLMT = MDEG + DELTA
             EHEAD = 0
C
  500        CONTINUE
                 MDNODE = DHEAD(MDEG)
                 IF  ( MDNODE .GT. 0 )  GO TO 600
                     MDEG = MDEG + 1
                     IF  ( MDEG .GT. MDLMT )  GO TO 900
                         GO TO 500
  600            CONTINUE
C                ----------------------------------------
C                REMOVE MDNODE FROM THE DEGREE STRUCTURE.
C                ----------------------------------------
                 NEXTMD = INVP(MDNODE)
                 DHEAD(MDEG) = NEXTMD
                 IF  ( NEXTMD .GT. 0 )  PERM(NEXTMD) = - MDEG
                 INVP(MDNODE) = - NUM
                 NOFSUB = NOFSUB + MDEG + QSIZE(MDNODE) - 2
                 IF  ( NUM+QSIZE(MDNODE) .GT. NEQNS )  GO TO 1000
C                ----------------------------------------------
C                ELIMINATE MDNODE AND PERFORM QUOTIENT GRAPH
C                TRANSFORMATION.  RESET TAG VALUE IF NECESSARY.
C                ----------------------------------------------
                 TAG = TAG + 1
                 IF  ( TAG .LT. MAXINT )  GO TO 800
                     TAG = 1
                     DO  700  I = 1, NEQNS
                         IF  ( MARKER(I) .LT. MAXINT )  MARKER(I) = 0
  700                CONTINUE
  800            CONTINUE
                 CALL  MMDELM ( MDNODE, XADJ, ADJNCY, DHEAD, INVP,
     1                          PERM, QSIZE, LLIST, MARKER, MAXINT,
     1                          TAG )
                 NUM = NUM + QSIZE(MDNODE)
                 LLIST(MDNODE) = EHEAD
                 EHEAD = MDNODE
                 IF  ( DELTA .GE. 0 )  GO TO 500
  900        CONTINUE
C            -------------------------------------------
C            UPDATE DEGREES OF THE NODES INVOLVED IN THE
C            MINIMUM DEGREE NODES ELIMINATION.
C            -------------------------------------------
             IF  ( NUM .GT. NEQNS )  GO TO 1000
             CALL  MMDUPD ( EHEAD, NEQNS, XADJ, ADJNCY, DELTA, MDEG,
     1                      DHEAD, INVP, PERM, QSIZE, LLIST, MARKER,
     1                      MAXINT, TAG )
             GO TO 300
C
 1000    CONTINUE
         CALL  MMDNUM ( NEQNS, PERM, INVP, QSIZE )
         RETURN
C
      END
C***************************************************************
C***************************************************************
C***     MMDINT ..... MULT MINIMUM DEGREE INITIALIZATION     ***
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE PERFORMS INITIALIZATION FOR THE
C        MULTIPLE ELIMINATION VERSION OF THE MINIMUM DEGREE
C        ALGORITHM.
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C        (XADJ,ADJNCY) - ADJACENCY STRUCTURE.
C
C     OUTPUT PARAMETERS -
C        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
C        QSIZE  - SIZE OF SUPERNODE (INITIALIZED TO ONE).
C        LLIST  - LINKED LIST.
C        MARKER - MARKER VECTOR.
C
C***************************************************************
C
      SUBROUTINE  MMDINT ( NEQNS, XADJ, ADJNCY, DHEAD, DFORW,
     1                     DBAKW, QSIZE, LLIST, MARKER )
C
C***************************************************************
C
C         INTEGER*2  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
C     1              LLIST(1) , MARKER(1), QSIZE(1)
         INTEGER*4  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
     1              LLIST(1) , MARKER(1), QSIZE(1)
         INTEGER*4  XADJ(1)
         INTEGER*4  FNODE , NDEG  , NEQNS , NODE
C
C***************************************************************
C
         DO  100  NODE = 1, NEQNS
             DHEAD(NODE) = 0
             QSIZE(NODE) = 1
             MARKER(NODE) = 0
             LLIST(NODE) = 0
  100    CONTINUE
C        ------------------------------------------
C        INITIALIZE THE DEGREE DOUBLY LINKED LISTS.
C        ------------------------------------------
         DO  200  NODE = 1, NEQNS
             NDEG = XADJ(NODE+1) - XADJ(NODE) + 1
             FNODE = DHEAD(NDEG)
             DFORW(NODE) = FNODE
             DHEAD(NDEG) = NODE
             IF  ( FNODE .GT. 0 )  DBAKW(FNODE) = NODE
             DBAKW(NODE) = - NDEG
  200    CONTINUE
         RETURN
C
      END
C***************************************************************
C***************************************************************
C**     MMDELM ..... MULTIPLE MINIMUM DEGREE ELIMINATION     ***
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE ELIMINATES THE NODE MDNODE OF
C        MINIMUM DEGREE FROM THE ADJACENCY STRUCTURE, WHICH
C        IS STORED IN THE QUOTIENT GRAPH FORMAT.  IT ALSO
C        TRANSFORMS THE QUOTIENT GRAPH REPRESENTATION OF THE
C        ELIMINATION GRAPH.
C
C     INPUT PARAMETERS -
C        MDNODE - NODE OF MINIMUM DEGREE.
C        MAXINT - ESTIMATE OF MAXIMUM REPRESENTABLE (SHORT)
C                 INTEGER.
C        TAG    - TAG VALUE.
C
C     UPDATED PARAMETERS -
C        (XADJ,ADJNCY) - UPDATED ADJACENCY STRUCTURE.
C        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
C        QSIZE  - SIZE OF SUPERNODE.
C        MARKER - MARKER VECTOR.
C        LLIST  - TEMPORARY LINKED LIST OF ELIMINATED NABORS.
C
C***************************************************************
C
      SUBROUTINE  MMDELM ( MDNODE, XADJ, ADJNCY, DHEAD, DFORW,
     1                     DBAKW, QSIZE, LLIST, MARKER, MAXINT,
     1                     TAG )
C
C***************************************************************
C
C         INTEGER*2  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
C     1              LLIST(1) , MARKER(1), QSIZE(1)
         INTEGER*4  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
     1              LLIST(1) , MARKER(1), QSIZE(1)
         INTEGER*4  XADJ(1)
         INTEGER*4  ELMNT , I     , ISTOP , ISTRT , J     ,
     1              JSTOP , JSTRT , LINK  , MAXINT, MDNODE,
     1              NABOR , NODE  , NPV   , NQNBRS, NXNODE,
     1              PVNODE, RLMT  , RLOC  , RNODE , TAG   ,
     1              XQNBR
C
C***************************************************************
C
C        -----------------------------------------------
C        FIND REACHABLE SET AND PLACE IN DATA STRUCTURE.
C        -----------------------------------------------
         MARKER(MDNODE) = TAG
         ISTRT = XADJ(MDNODE)
         ISTOP = XADJ(MDNODE+1) - 1
C        -------------------------------------------------------
C        ELMNT POINTS TO THE BEGINNING OF THE LIST OF ELIMINATED
C        NABORS OF MDNODE, AND RLOC GIVES THE STORAGE LOCATION
C        FOR THE NEXT REACHABLE NODE.
C        -------------------------------------------------------
         ELMNT = 0
         RLOC = ISTRT
         RLMT = ISTOP
         DO  200  I = ISTRT, ISTOP
             NABOR = ADJNCY(I)
             IF  ( NABOR .EQ. 0 )  GO TO 300
                 IF  ( MARKER(NABOR) .GE. TAG )  GO TO 200
                     MARKER(NABOR) = TAG
                     IF  ( DFORW(NABOR) .LT. 0 )  GO TO 100
                         ADJNCY(RLOC) = NABOR
                         RLOC = RLOC + 1
                         GO TO 200
  100                CONTINUE
                     LLIST(NABOR) = ELMNT
                     ELMNT = NABOR
  200    CONTINUE
  300    CONTINUE
C            -----------------------------------------------------
C            MERGE WITH REACHABLE NODES FROM GENERALIZED ELEMENTS.
C            -----------------------------------------------------
             IF  ( ELMNT .LE. 0 )  GO TO 1000
                 ADJNCY(RLMT) = - ELMNT
                 LINK = ELMNT
  400            CONTINUE
                     JSTRT = XADJ(LINK)
                     JSTOP = XADJ(LINK+1) - 1
                     DO  800  J = JSTRT, JSTOP
                         NODE = ADJNCY(J)
                         LINK = - NODE
                         IF  ( NODE )  400, 900, 500
  500                    CONTINUE
                         IF  ( MARKER(NODE) .GE. TAG  .OR.
     1                         DFORW(NODE) .LT. 0 )  GO TO 800
                             MARKER(NODE) = TAG
C                            ---------------------------------
C                            USE STORAGE FROM ELIMINATED NODES
C                            IF NECESSARY.
C                            ---------------------------------
  600                        CONTINUE
                                 IF  ( RLOC .LT. RLMT )  GO TO 700
                                     LINK = - ADJNCY(RLMT)
                                     RLOC = XADJ(LINK)
                                     RLMT = XADJ(LINK+1) - 1
                                     GO TO 600
  700                        CONTINUE
                             ADJNCY(RLOC) = NODE
                             RLOC = RLOC + 1
  800                CONTINUE
  900            CONTINUE
                 ELMNT = LLIST(ELMNT)
                 GO TO 300
 1000    CONTINUE
         IF  ( RLOC .LE. RLMT )  ADJNCY(RLOC) = 0
C        --------------------------------------------------------
C        FOR EACH NODE IN THE REACHABLE SET, DO THE FOLLOWING ...
C        --------------------------------------------------------
         LINK = MDNODE
 1100    CONTINUE
             ISTRT = XADJ(LINK)
             ISTOP = XADJ(LINK+1) - 1
             DO  1700  I = ISTRT, ISTOP
                 RNODE = ADJNCY(I)
                 LINK = - RNODE
                 IF  ( RNODE )  1100, 1800, 1200
 1200            CONTINUE
C                --------------------------------------------
C                IF RNODE IS IN THE DEGREE LIST STRUCTURE ...
C                --------------------------------------------
                 PVNODE = DBAKW(RNODE)
                 IF  ( PVNODE .EQ. 0  .OR.
     1                 PVNODE .EQ. (-MAXINT) )  GO TO 1300
C                    -------------------------------------
C                    THEN REMOVE RNODE FROM THE STRUCTURE.
C                    -------------------------------------
                     NXNODE = DFORW(RNODE)
                     IF  ( NXNODE .GT. 0 )  DBAKW(NXNODE) = PVNODE
                     IF  ( PVNODE .GT. 0 )  DFORW(PVNODE) = NXNODE
                     NPV = - PVNODE
                     IF  ( PVNODE .LT. 0 )  DHEAD(NPV) = NXNODE
 1300            CONTINUE
C                ----------------------------------------
C                PURGE INACTIVE QUOTIENT NABORS OF RNODE.
C                ----------------------------------------
                 JSTRT = XADJ(RNODE)
                 JSTOP = XADJ(RNODE+1) - 1
                 XQNBR = JSTRT
                 DO  1400  J = JSTRT, JSTOP
                     NABOR = ADJNCY(J)
                     IF  ( NABOR .EQ. 0 )  GO TO 1500
                         IF  ( MARKER(NABOR) .GE. TAG )  GO TO 1400
                             ADJNCY(XQNBR) = NABOR
                             XQNBR = XQNBR + 1
 1400            CONTINUE
 1500            CONTINUE
C                ----------------------------------------
C                IF NO ACTIVE NABOR AFTER THE PURGING ...
C                ----------------------------------------
                 NQNBRS = XQNBR - JSTRT
                 IF  ( NQNBRS .GT. 0 )  GO TO 1600
C                    -----------------------------
C                    THEN MERGE RNODE WITH MDNODE.
C                    -----------------------------
                     QSIZE(MDNODE) = QSIZE(MDNODE) + QSIZE(RNODE)
                     QSIZE(RNODE) = 0
                     MARKER(RNODE) = MAXINT
                     DFORW(RNODE) = - MDNODE
                     DBAKW(RNODE) = - MAXINT
                     GO TO 1700
 1600            CONTINUE
C                --------------------------------------
C                ELSE FLAG RNODE FOR DEGREE UPDATE, AND
C                ADD MDNODE AS A NABOR OF RNODE.
C                --------------------------------------
                 DFORW(RNODE) = NQNBRS + 1
                 DBAKW(RNODE) = 0
                 ADJNCY(XQNBR) = MDNODE
                 XQNBR = XQNBR + 1
                 IF  ( XQNBR .LE. JSTOP )  ADJNCY(XQNBR) = 0
C
 1700        CONTINUE
 1800    CONTINUE
         RETURN
C
      END
C***************************************************************
C***************************************************************
C*****     MMDUPD ..... MULTIPLE MINIMUM DEGREE UPDATE     *****
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE UPDATES THE DEGREES OF NODES
C        AFTER A MULTIPLE ELIMINATION STEP.
C
C     INPUT PARAMETERS -
C        EHEAD  - THE BEGINNING OF THE LIST OF ELIMINATED
C                 NODES (I.E., NEWLY FORMED ELEMENTS).
C        NEQNS  - NUMBER OF EQUATIONS.
C        (XADJ,ADJNCY) - ADJACENCY STRUCTURE.
C        DELTA  - TOLERANCE VALUE FOR MULTIPLE ELIMINATION.
C        MAXINT - MAXIMUM MACHINE REPRESENTABLE (SHORT)
C                 INTEGER.
C
C     UPDATED PARAMETERS -
C        MDEG   - NEW MINIMUM DEGREE AFTER DEGREE UPDATE.
C        (DHEAD,DFORW,DBAKW) - DEGREE DOUBLY LINKED STRUCTURE.
C        QSIZE  - SIZE OF SUPERNODE.
C        LLIST  - WORKING LINKED LIST.
C        MARKER - MARKER VECTOR FOR DEGREE UPDATE.
C        TAG    - TAG VALUE.
C
C***************************************************************
C
      SUBROUTINE  MMDUPD ( EHEAD, NEQNS, XADJ, ADJNCY, DELTA,
     1                     MDEG, DHEAD, DFORW, DBAKW, QSIZE,
     1                     LLIST, MARKER, MAXINT, TAG )
C
C***************************************************************
C
C         INTEGER*2  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
C     1              LLIST(1) , MARKER(1), QSIZE(1)
         INTEGER*4  ADJNCY(1), DBAKW(1) , DFORW(1) , DHEAD(1) ,
     1              LLIST(1) , MARKER(1), QSIZE(1)
         INTEGER*4  XADJ(1)
         INTEGER*4  DEG   , DEG0  , DELTA , EHEAD , ELMNT ,
     1              ENODE , FNODE , I     , IQ2   , ISTOP ,
     1              ISTRT , J     , JSTOP , JSTRT , LINK  ,
     1              MAXINT, MDEG  , MDEG0 , MTAG  , NABOR ,
     1              NEQNS , NODE  , Q2HEAD, QXHEAD, TAG
C
C***************************************************************
C
         MDEG0 = MDEG + DELTA
         ELMNT = EHEAD
  100    CONTINUE
C            -------------------------------------------------------
C            FOR EACH OF THE NEWLY FORMED ELEMENT, DO THE FOLLOWING.
C            (RESET TAG VALUE IF NECESSARY.)
C            -------------------------------------------------------
             IF  ( ELMNT .LE. 0 )  RETURN
             MTAG = TAG + MDEG0
             IF  ( MTAG .LT. MAXINT )  GO TO 300
                 TAG = 1
                 DO  200  I = 1, NEQNS
                     IF  ( MARKER(I) .LT. MAXINT )  MARKER(I) = 0
  200            CONTINUE
                 MTAG = TAG + MDEG0
  300        CONTINUE
C            ---------------------------------------------
C            CREATE TWO LINKED LISTS FROM NODES ASSOCIATED
C            WITH ELMNT: ONE WITH TWO NABORS (Q2HEAD) IN
C            ADJACENCY STRUCTURE, AND THE OTHER WITH MORE
C            THAN TWO NABORS (QXHEAD).  ALSO COMPUTE DEG0,
C            NUMBER OF NODES IN THIS ELEMENT.
C            ---------------------------------------------
             Q2HEAD = 0
             QXHEAD = 0
             DEG0 = 0
             LINK = ELMNT
  400        CONTINUE
                 ISTRT = XADJ(LINK)
                 ISTOP = XADJ(LINK+1) - 1
                 DO  700  I = ISTRT, ISTOP
                     ENODE = ADJNCY(I)
                     LINK = - ENODE
                     IF  ( ENODE )  400, 800, 500
C
  500                CONTINUE
                     IF  ( QSIZE(ENODE) .EQ. 0 )  GO TO 700
                         DEG0 = DEG0 + QSIZE(ENODE)
                         MARKER(ENODE) = MTAG
C                        ----------------------------------
C                        IF ENODE REQUIRES A DEGREE UPDATE,
C                        THEN DO THE FOLLOWING.
C                        ----------------------------------
                         IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 700
C                            ---------------------------------------
C                            PLACE EITHER IN QXHEAD OR Q2HEAD LISTS.
C                            ---------------------------------------
                             IF  ( DFORW(ENODE) .EQ. 2 )  GO TO 600
                                 LLIST(ENODE) = QXHEAD
                                 QXHEAD = ENODE
                                 GO TO 700
  600                        CONTINUE
                             LLIST(ENODE) = Q2HEAD
                             Q2HEAD = ENODE
  700            CONTINUE
  800        CONTINUE
C            --------------------------------------------
C            FOR EACH ENODE IN Q2 LIST, DO THE FOLLOWING.
C            --------------------------------------------
             ENODE = Q2HEAD
             IQ2 = 1
  900        CONTINUE
                 IF  ( ENODE .LE. 0 )  GO TO 1500
                 IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 2200
                     TAG = TAG + 1
                     DEG = DEG0
C                    ------------------------------------------
C                    IDENTIFY THE OTHER ADJACENT ELEMENT NABOR.
C                    ------------------------------------------
                     ISTRT = XADJ(ENODE)
                     NABOR = ADJNCY(ISTRT)
                     IF  ( NABOR .EQ. ELMNT )  NABOR = ADJNCY(ISTRT+1)
C                    ------------------------------------------------
C                    IF NABOR IS UNELIMINATED, INCREASE DEGREE COUNT.
C                    ------------------------------------------------
                     LINK = NABOR
                     IF  ( DFORW(NABOR) .LT. 0 )  GO TO 1000
                         DEG = DEG + QSIZE(NABOR)
                         GO TO 2100
 1000                CONTINUE
C                        --------------------------------------------
C                        OTHERWISE, FOR EACH NODE IN THE 2ND ELEMENT,
C                        DO THE FOLLOWING.
C                        --------------------------------------------
                         ISTRT = XADJ(LINK)
                         ISTOP = XADJ(LINK+1) - 1
                         DO  1400  I = ISTRT, ISTOP
                             NODE = ADJNCY(I)
                             LINK = - NODE
                             IF  ( NODE .EQ. ENODE )  GO TO 1400
                             IF  ( NODE )  1000, 2100, 1100
C
 1100                        CONTINUE
                             IF  ( QSIZE(NODE) .EQ. 0 )  GO TO 1400
                             IF  ( MARKER(NODE) .GE. TAG )  GO TO 1200
C                                -------------------------------------
C                                CASE WHEN NODE IS NOT YET CONSIDERED.
C                                -------------------------------------
                                 MARKER(NODE) = TAG
                                 DEG = DEG + QSIZE(NODE)
                                 GO TO 1400
 1200                        CONTINUE
C                            ----------------------------------------
C                            CASE WHEN NODE IS INDISTINGUISHABLE FROM
C                            ENODE.  MERGE THEM INTO A NEW SUPERNODE.
C                            ----------------------------------------
                             IF  ( DBAKW(NODE) .NE. 0 )  GO TO 1400
                             IF  ( DFORW(NODE) .NE. 2 )  GO TO 1300
                                 QSIZE(ENODE) = QSIZE(ENODE) +
     1                                          QSIZE(NODE)
                                 QSIZE(NODE) = 0
                                 MARKER(NODE) = MAXINT
                                 DFORW(NODE) = - ENODE
                                 DBAKW(NODE) = - MAXINT
                                 GO TO 1400
 1300                        CONTINUE
C                            --------------------------------------
C                            CASE WHEN NODE IS OUTMATCHED BY ENODE.
C                            --------------------------------------
                             IF  ( DBAKW(NODE) .EQ.0 )
     1                             DBAKW(NODE) = - MAXINT
 1400                    CONTINUE
                         GO TO 2100
 1500            CONTINUE
C                ------------------------------------------------
C                FOR EACH ENODE IN THE QX LIST, DO THE FOLLOWING.
C                ------------------------------------------------
                 ENODE = QXHEAD
                 IQ2 = 0
 1600            CONTINUE
                     IF  ( ENODE .LE. 0 )  GO TO 2300
                     IF  ( DBAKW(ENODE) .NE. 0 )  GO TO 2200
                         TAG = TAG + 1
                         DEG = DEG0
C                        ---------------------------------
C                        FOR EACH UNMARKED NABOR OF ENODE,
C                        DO THE FOLLOWING.
C                        ---------------------------------
                         ISTRT = XADJ(ENODE)
                         ISTOP = XADJ(ENODE+1) - 1
                         DO  2000  I = ISTRT, ISTOP
                             NABOR = ADJNCY(I)
                             IF  ( NABOR .EQ. 0 )  GO TO 2100
                             IF  ( MARKER(NABOR) .GE. TAG )  GO TO 2000
                                 MARKER(NABOR) = TAG
                                 LINK = NABOR
C                                ------------------------------
C                                IF UNELIMINATED, INCLUDE IT IN
C                                DEG COUNT.
C                                ------------------------------
                                 IF  ( DFORW(NABOR) .LT. 0 )  GO TO 1700
                                     DEG = DEG + QSIZE(NABOR)
                                     GO TO 2000
 1700                            CONTINUE
C                                    -------------------------------
C                                    IF ELIMINATED, INCLUDE UNMARKED
C                                    NODES IN THIS ELEMENT INTO THE
C                                    DEGREE COUNT.
C                                    -------------------------------
                                     JSTRT = XADJ(LINK)
                                     JSTOP = XADJ(LINK+1) - 1
                                     DO  1900  J = JSTRT, JSTOP
                                         NODE = ADJNCY(J)
                                         LINK = - NODE
                                         IF  ( NODE )  1700, 2000, 1800
C
 1800                                    CONTINUE
                                         IF  ( MARKER(NODE) .GE. TAG )
     1                                         GO TO 1900
                                             MARKER(NODE) = TAG
                                             DEG = DEG + QSIZE(NODE)
 1900                                CONTINUE
 2000                    CONTINUE
 2100                CONTINUE
C                    -------------------------------------------
C                    UPDATE EXTERNAL DEGREE OF ENODE IN DEGREE
C                    STRUCTURE, AND MDEG (MIN DEG) IF NECESSARY.
C                    -------------------------------------------
                     DEG = DEG - QSIZE(ENODE) + 1
                     FNODE = DHEAD(DEG)
                     DFORW(ENODE) = FNODE
                     DBAKW(ENODE) = - DEG
                     IF  ( FNODE .GT. 0 )  DBAKW(FNODE) = ENODE
                     DHEAD(DEG) = ENODE
                     IF  ( DEG .LT. MDEG )  MDEG = DEG
 2200                CONTINUE
C                    ----------------------------------
C                    GET NEXT ENODE IN CURRENT ELEMENT.
C                    ----------------------------------
                     ENODE = LLIST(ENODE)
                     IF  ( IQ2 .EQ. 1 )  GO TO 900
                         GO TO 1600
 2300        CONTINUE
C            -----------------------------
C            GET NEXT ELEMENT IN THE LIST.
C            -----------------------------
             TAG = MTAG
             ELMNT = LLIST(ELMNT)
             GO TO 100
C
      END
C***************************************************************
C***************************************************************
C*****     MMDNUM ..... MULTI MINIMUM DEGREE NUMBERING     *****
C***************************************************************
C***************************************************************
C
C     AUTHOR - JOSEPH W.H. LIU
C              DEPT OF COMPUTER SCIENCE, YORK UNIVERSITY.
C
C     PURPOSE - THIS ROUTINE PERFORMS THE FINAL STEP IN
C        PRODUCING THE PERMUTATION AND INVERSE PERMUTATION
C        VECTORS IN THE MULTIPLE ELIMINATION VERSION OF THE
C        MINIMUM DEGREE ORDERING ALGORITHM.
C
C     INPUT PARAMETERS -
C        NEQNS  - NUMBER OF EQUATIONS.
C        QSIZE  - SIZE OF SUPERNODES AT ELIMINATION.
C
C     UPDATED PARAMETERS -
C        INVP   - INVERSE PERMUTATION VECTOR.  ON INPUT,
C                 IF QSIZE(NODE)=0, THEN NODE HAS BEEN MERGED
C                 INTO THE NODE -INVP(NODE); OTHERWISE,
C                 -INVP(NODE) IS ITS INVERSE LABELLING.
C
C     OUTPUT PARAMETERS -
C        PERM   - THE PERMUTATION VECTOR.
C
C***************************************************************
C
      SUBROUTINE  MMDNUM ( NEQNS, PERM, INVP, QSIZE )
C
C***************************************************************
C
C         INTEGER*2  INVP(1)  , PERM(1)  , QSIZE(1)
         INTEGER*4  INVP(1)  , PERM(1)  , QSIZE(1)
         INTEGER*4  FATHER, NEQNS , NEXTF , NODE  , NQSIZE,
     1              NUM   , ROOT
C
C***************************************************************
C
         DO  100  NODE = 1, NEQNS
             NQSIZE = QSIZE(NODE)
             IF  ( NQSIZE .LE. 0 )  PERM(NODE) = INVP(NODE)
             IF  ( NQSIZE .GT. 0 )  PERM(NODE) = - INVP(NODE)
  100    CONTINUE
C        ------------------------------------------------------
C        FOR EACH NODE WHICH HAS BEEN MERGED, DO THE FOLLOWING.
C        ------------------------------------------------------
         DO  500  NODE = 1, NEQNS
             IF  ( PERM(NODE) .GT. 0 )  GO TO 500
C                -----------------------------------------
C                TRACE THE MERGED TREE UNTIL ONE WHICH HAS
C                NOT BEEN MERGED, CALL IT ROOT.
C                -----------------------------------------
                 FATHER = NODE
  200            CONTINUE
                     IF  ( PERM(FATHER) .GT. 0 )  GO TO 300
                         FATHER = - PERM(FATHER)
                         GO TO 200
  300            CONTINUE
C                -----------------------
C                NUMBER NODE AFTER ROOT.
C                -----------------------
                 ROOT = FATHER
                 NUM = PERM(ROOT) + 1
                 INVP(NODE) = - NUM
                 PERM(ROOT) = NUM
C                ------------------------
C                SHORTEN THE MERGED TREE.
C                ------------------------
                 FATHER = NODE
  400            CONTINUE
                     NEXTF = - PERM(FATHER)
                     IF  ( NEXTF .LE. 0 )  GO TO 500
                         PERM(FATHER) = - ROOT
                         FATHER = NEXTF
                         GO TO 400
  500    CONTINUE
C        ----------------------
C        READY TO COMPUTE PERM.
C        ----------------------
         DO  600  NODE = 1, NEQNS
             NUM = - INVP(NODE)
             INVP(NODE) = NUM
             PERM(NUM) = NODE
  600    CONTINUE
         RETURN
C
      END


