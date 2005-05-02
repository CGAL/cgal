c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

        subroutine ireadhb(fname,type,
     $                     nrows,ncols,nnz)

      
        implicit none
        integer nrows, ncols, nnz
        integer totcrd, ptrcrd,
     $		indcrd, valcrd, rhscrd
c        character title*72, key*30, type*3, ptrfmt*16,
c     $          indfmt*16, valfmt*20, rhsfmt*20
        character title*72, key*30, type*3
        character fname*256
c        logical sym
c        double precision skew
c        double precision myrand
c        character rhstyp*3
c        integer nzrhs
        integer nel

c-----------------------------------------------------------------------

c       read header information from Harwell/Boeing matrix
        open (99, file=fname, err=999, status="OLD")

        read (99, 10, err = 999)
     $          title, key,
     $          totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     $          type, nrows, ncols, nnz, nel
10      format (a72, a8 / 5i14 / a3, 11x, 4i14)

        write (0, 30)
     $          title, key,
     $          totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     $          type, nrows, ncols, nnz, nel
30      format (
     $          ' title: ', a72 /
     $          ' key: ', a8 /
     $          ' Lines: tot: ', i14,' ptr: ',i14,' ind: ',i14 /
     $          '        val: ', i14,' rhs: ',i14 /
     $          ' type: ', a3, ' nrow: ', i14, ' ncol: ', i14 /
     $          ' nz: ', i14, ' elements: ', i14)

        close (99)
        return

999     write (0,*) 'Read error: Harwell/Boeing matrix'
        stop
        end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

        subroutine dreadhb(fname,
     $                     nrows,ncols,nnz,
     $                     Ptr,Index,Value)

        implicit none
        integer nrows, ncols, nnz
        integer Ptr(*), Index(*), totcrd, ptrcrd,
     $		indcrd, valcrd, rhscrd, nrhs, row, col, p
        character title*72, key*30, type*3, ptrfmt*16,
     $          indfmt*16, valfmt*20, rhsfmt*20
        character fname*256
        logical sym
        double precision Value (*), skew
        double precision myrand
        character rhstyp*3
        integer nzrhs, nel

c-----------------------------------------------------------------------

c       read header information from Harwell/Boeing matrix

        open (99, file=fname, err=198, status="OLD")

        read (99, 105, err = 198)
     $          title, key,
     $          totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     $          type, nrows, ncols, nnz, nel
105     format (a72, a8 / 5i14 / a3, 11x, 4i14)

        read (99, 110, err = 198)
     $          ptrfmt, indfmt, valfmt, rhsfmt
        if (rhscrd .gt. 0) then
c          new Harwell/Boeing format:
           read (99, 120, err = 198) rhstyp,nrhs,nzrhs
           endif
110     format (2a16, 2a20)
120     format (a3, 11x, 2i14)

        skew = 0.0
        if (type (2:2) .eq. 'Z' .or. type (2:2) .eq. 'z') skew = -1.0
        if (type (2:2) .eq. 'S' .or. type (2:2) .eq. 's') skew =  1.0
        sym = skew .ne. 0.0

        write (0, 130)
     $          ptrfmt, indfmt, valfmt, rhsfmt
        if (rhscrd .gt. 0) then
c          new Harwell/Boeing format:
           write (0, 140) rhstyp,nrhs,nzrhs
           endif
130     format (
     $          ' ptrfmt: ', a20, ' rowfmt: ', a20, /
     $          ' valfmt: ', a20, ' rhsfmt: ', a20)
140     format (' rhstyp: ', a3, ' nrhs: ', i14, ' nzrhs: ', i14)
        write (0, *) ' sym: ', sym, ' skew: ', skew

        print *,'reading colptr'

        read (99, ptrfmt, err = 198) (Ptr (p), p = 1, ncols+1)
        print *,'reading rowind'
        read (99, indfmt, err = 198) (Index (p), p = 1, nnz)

c      what's this? maybe for rectangualr matrices

c        do 155 col = ncols+2, ncols+1
c           Ptr (col) = Ptr (ncols+1)
c155         continue

        print *,'reading values'
c       read the values, or create random-valued matrix
        if (valcrd .gt. 0) then
           read (99, valfmt, err = 198) (Value (p), p = 1, nnz)
        else
          if (sym) then
            do 157 col = 1, ncols
              do 156 p = Ptr(col), Ptr(col+1)-1
                row = Index(p)
                if (row .eq. col) then
                  Value(p) = ncols
                else
                  Value(p) = -1.0
                endif
156           continue
157         continue
          else
            Value (1) = myrand (0)
            do 158 p = 1, nnz
               Value (p) = myrand (-1)
158         continue
          endif
        endif

c  create the triplet form of the input matrix

c        do 100 col = 1, n
c           do 90 p = Ptr (col), Ptr (col+1) - 1
c              row = Index (p)
c              write (6, 200) row, col, Value (p)
c              if (sym .and. row .ne. col) then
c		 write (6, 200) col, row, skew * Value (p)
c		 endif
c90            continue
c100        continue
c200	format (2i7, e26.16e3)

        close (99)
        return

198     write (0,*) 'Read error: Harwell/Boeing matrix'
        stop
        end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

        subroutine sreadhb(fname,
     $                     nrows,ncols,nnz,
     $                     Ptr,Index,Value)

        implicit none
        integer nrows, ncols, nnz
        integer Ptr(*), Index(*), totcrd, ptrcrd,
     $		indcrd, valcrd, rhscrd, nrhs, row, col, p
        character title*72, key*30, type*3, ptrfmt*16,
     $          indfmt*16, valfmt*20, rhsfmt*20
        character fname*256
        logical sym
        real*4 Value (*), skew
        double precision myrand
        character rhstyp*3
        integer nzrhs, nel

c-----------------------------------------------------------------------

c       read header information from Harwell/Boeing matrix

        open (99, file=fname, err=298, status="OLD")

        read (99, 205, err = 298)
     $          title, key,
     $          totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     $          type, nrows, ncols, nnz, nel
205     format (a72, a8 / 5i14 / a3, 11x, 4i14)

        read (99, 210, err = 298)
     $          ptrfmt, indfmt, valfmt, rhsfmt
        if (rhscrd .gt. 0) then
c          new Harwell/Boeing format:
           read (99, 220, err = 298) rhstyp,nrhs,nzrhs
           endif
210     format (2a16, 2a20)
220     format (a3, 11x, 2i14)

        skew = 0.0
        if (type (2:2) .eq. 'Z' .or. type (2:2) .eq. 'z') skew = -1.0
        if (type (2:2) .eq. 'S' .or. type (2:2) .eq. 's') skew =  1.0
        sym = skew .ne. 0.0

        write (0, 230)
     $          ptrfmt, indfmt, valfmt, rhsfmt
        if (rhscrd .gt. 0) then
c          new Harwell/Boeing format:
           write (0, 240) rhstyp,nrhs,nzrhs
           endif
230     format (
     $          ' ptrfmt: ', a20, ' rowfmt: ', a20, /
     $          ' valfmt: ', a20, ' rhsfmt: ', a20)
240     format (' rhstyp: ', a3, ' nrhs: ', i14, ' nzrhs: ', i14)
        write (0, *) ' sym: ', sym, ' skew: ', skew

        print *,'reading colptr'

        read (99, ptrfmt, err = 298) (Ptr (p), p = 1, ncols+1)
        print *,'reading rowind'
        read (99, indfmt, err = 298) (Index (p), p = 1, nnz)

c      what's this? maybe for rectangualr matrices

c        do 255 col = ncols+2, ncols+1
c           Ptr (col) = Ptr (ncols+1)
c255         continue

        print *,'reading values'
c       read the values, or create random-valued matrix
        if (valcrd .gt. 0) then
           read (99, valfmt, err = 298) (Value (p), p = 1, nnz)
        else
          if (sym) then
            do 257 col = 1, ncols
              do 256 p = Ptr(col), Ptr(col+1)-1
                row = Index(p)
                if (row .eq. col) then
                  Value(p) = ncols
                else
                  Value(p) = -1.0
                endif
256           continue
257         continue
          else
            Value (1) = myrand (0)
            do 258 p = 1, nnz
               Value (p) = myrand (-1)
258         continue
          endif
        endif

c  create the triplet form of the input matrix

c        do 100 col = 1, n
c           do 90 p = Ptr (col), Ptr (col+1) - 1
c              row = Index (p)
c              write (6, 200) row, col, Value (p)
c              if (sym .and. row .ne. col) then
c		 write (6, 200) col, row, skew * Value (p)
c		 endif
c90            continue
c100        continue
c200	format (2i7, e26.16e3)

        close (99)
        return

298     write (0,*) 'Read error: Harwell/Boeing matrix'
        stop
        end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

        subroutine zreadhb(fname,
     $                     nrows,ncols,nnz,
     $                     Ptr,Index,Value)

        implicit none
        integer nrows, ncols, nnz
        integer Ptr(*), Index(*), totcrd, ptrcrd,
     $		indcrd, valcrd, rhscrd, nrhs, row, col, p
        character title*72, key*30, type*3, ptrfmt*16,
     $          indfmt*16, valfmt*20, rhsfmt*20
        character fname*256
        logical sym
        double complex Value (*), skew
        double precision myrand
        character rhstyp*3
        integer nzrhs, nel

c-----------------------------------------------------------------------

c       read header information from Harwell/Boeing matrix

        open (99, file=fname, err=398, status="OLD")

        read (99, 305, err = 398)
     $          title, key,
     $          totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     $          type, nrows, ncols, nnz, nel
305     format (a72, a8 / 5i14 / a3, 11x, 4i14)

        read (99, 310, err = 398)
     $          ptrfmt, indfmt, valfmt, rhsfmt
        if (rhscrd .gt. 0) then
c          new Harwell/Boeing format:
           read (99, 320, err = 398) rhstyp,nrhs,nzrhs
           endif
310      format (2a16, 2a20)
320      format (a3, 11x, 2i14)

        skew = 0.0
        if (type (2:2) .eq. 'Z' .or. type (2:2) .eq. 'z') skew = -1.0
        if (type (2:2) .eq. 'S' .or. type (2:2) .eq. 's') skew =  1.0
        sym = skew .ne. 0.0

        write (0, 330)
     $          ptrfmt, indfmt, valfmt, rhsfmt
        if (rhscrd .gt. 0) then
c          new Harwell/Boeing format:
           write (0, 340) rhstyp,nrhs,nzrhs
           endif
330     format (
     $          ' ptrfmt: ', a20, ' rowfmt: ', a20, /
     $          ' valfmt: ', a20, ' rhsfmt: ', a20)
340     format (' rhstyp: ', a3, ' nrhs: ', i14, ' nzrhs: ', i14)
        write (0, *) ' sym: ', sym, ' skew: ', skew

        print *,'reading colptr'

        read (99, ptrfmt, err = 398) (Ptr (p), p = 1, ncols+1)
        print *,'reading rowind'
        read (99, indfmt, err = 398) (Index (p), p = 1, nnz)

c      what's this? maybe for rectangualr matrices

c        do 355 col = ncols+2, ncols+1
c           Ptr (col) = Ptr (ncols+1)
c355        continue

        print *,'reading values'
c       read the values, or create random-valued matrix
        if (valcrd .gt. 0) then
           read (99, valfmt, err = 398) (Value (p), p = 1, nnz)
        else
          if (sym) then
            do 357 col = 1, ncols
              do 356 p = Ptr(col), Ptr(col+1)-1
                row = Index(p)
                if (row .eq. col) then
                  Value(p) = ncols
                else
                  Value(p) = -1.0
                endif
356           continue
357         continue
          else
            Value (1) = myrand (0)
            do 350 p = 1, nnz
               Value (p) = myrand (-1)
350         continue
          endif
        endif

c  create the triplet form of the input matrix

c        do 100 col = 1, n
c           do 90 p = Ptr (col), Ptr (col+1) - 1
c              row = Index (p)
c              write (6, 200) row, col, Value (p)
c              if (sym .and. row .ne. col) then
c		 write (6, 200) col, row, skew * Value (p)
c		 endif
c90            continue
c100        continue
c200	format (2i7, e26.16e3)

        close (99)
        return

398     write (0,*) 'Read error: Harwell/Boeing matrix'
        stop
        end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

        subroutine creadhb(fname,
     $                     nrows,ncols,nnz,
     $                     Ptr,Index,Value)

        implicit none
        integer nrows, ncols, nnz
        integer Ptr(*), Index(*), totcrd, ptrcrd,
     $		indcrd, valcrd, rhscrd, nrhs, row, col, p
        character title*72, key*30, type*3, ptrfmt*16,
     $          indfmt*16, valfmt*20, rhsfmt*20
        character fname*256
        logical sym
        complex Value (*), skew
        double precision myrand
        character rhstyp*3
        integer nzrhs, nel

c-----------------------------------------------------------------------

c       read header information from Harwell/Boeing matrix

        open (99, file=fname, err=498, status="OLD")

        read (99, 405, err = 498)
     $          title, key,
     $          totcrd, ptrcrd, indcrd, valcrd, rhscrd,
     $          type, nrows, ncols, nnz, nel
405     format (a72, a8 / 5i14 / a3, 11x, 4i14)

        read (99, 410, err = 498)
     $          ptrfmt, indfmt, valfmt, rhsfmt
        if (rhscrd .gt. 0) then
c          new Harwell/Boeing format:
           read (99, 420, err = 498) rhstyp,nrhs,nzrhs
           endif
410     format (2a16, 2a20)
420     format (a3, 11x, 2i14)

        skew = 0.0
        if (type (2:2) .eq. 'Z' .or. type (2:2) .eq. 'z') skew = -1.0
        if (type (2:2) .eq. 'S' .or. type (2:2) .eq. 's') skew =  1.0
        sym = skew .ne. 0.0

        write (0, 430)
     $          ptrfmt, indfmt, valfmt, rhsfmt
        if (rhscrd .gt. 0) then
c          new Harwell/Boeing format:
           write (0, 440) rhstyp,nrhs,nzrhs
           endif
430     format (
     $          ' ptrfmt: ', a20, ' rowfmt: ', a20, /
     $          ' valfmt: ', a20, ' rhsfmt: ', a20)
440     format (' rhstyp: ', a3, ' nrhs: ', i14, ' nzrhs: ', i14)
        write (0, *) ' sym: ', sym, ' skew: ', skew

        print *,'reading colptr'

        read (99, ptrfmt, err = 498) (Ptr (p), p = 1, ncols+1)
        print *,'reading rowind'
        read (99, indfmt, err = 498) (Index (p), p = 1, nnz)

c      what's this? maybe for rectangualr matrices

c        do 455 col = ncols+2, ncols+1
c           Ptr (col) = Ptr (ncols+1)
c455        continue

        print *,'reading values'
c       read the values, or create random-valued matrix
        if (valcrd .gt. 0) then
           read (99, valfmt, err = 498) (Value (p), p = 1, nnz)
        else
          if (sym) then
            do 457 col = 1, ncols
              do 456 p = Ptr(col), Ptr(col+1)-1
                row = Index(p)
                if (row .eq. col) then
                  Value(p) = ncols
                else
                  Value(p) = -1.0
                endif
456           continue
457         continue
          else
            Value (1) = myrand (0)
            do 450 p = 1, nnz
               Value (p) = myrand (-1)
450          continue
          endif
        endif

c  create the triplet form of the input matrix

c        do 100 col = 1, n
c           do 90 p = Ptr (col), Ptr (col+1) - 1
c              row = Index (p)
c              write (6, 200) row, col, Value (p)
c              if (sym .and. row .ne. col) then
c		 write (6, 200) col, row, skew * Value (p)
c		 endif
c90            continue
c100        continue
c200	format (2i7, e26.16e3)

        close (99)
        return

498     write (0,*) 'Read error: Harwell/Boeing matrix'
        stop
        end

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

c=== Myrand ============================================================
c
c  Derived from the FA01 routines in the MUPS package (CERFACS and/or
c  Harwell).  CERFACS and/or Harwell copyrights may apply.  Permission
c  granted to use this routine in the DEMO PROGRAM only.
c
c  DEMO PROGRAM.
c
c  random number generator
c  i = 0:  reinitialize the sequence
c  i >=0:  return 0 < x < 1
c  i < 0:  return -1 < x < 1

        double precision function myrand (i)
        integer i
        double precision seed, start, mfac, d2to32
        common /mrand/ seed
        parameter (start = 1431655765.d0,
     $             d2to32 = 4294967296.d0, mfac = 9228907.d0)

        if (i .eq. 0) then
c          reinitialize to known sequence
           seed = start
           endif
        seed = dmod (seed * mfac, d2to32)

        if (i .ge. 0) then
           myrand = (seed/d2to32)
        else
           myrand = 2 * (seed/d2to32) - 1
           endif
        return
        end




