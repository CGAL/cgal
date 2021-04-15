// Copyright (c) 2005-2008 ASCLEPIOS Project, INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of the ImageIO Library, and as been adapted for CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     :  ASCLEPIOS Project (INRIA Sophia-Antipolis), Laurent Rineau

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#ifdef _MSC_VER
#pragma warning ( disable : 4068 4786 4081 )
#endif

#ifdef MINC_FILES

#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>



/** Magic header for MINC (MNI NetCDF) file format */
#ifdef MINC_FILES
#  define MINC_MAGIC "CDF"
#endif


CGAL_INLINE_FUNCTION
int readMincHeader(_image *im, const char* filename,
                   double *startx, double *starty, double *startz,
                   double *stepx, double *stepy, double *stepz,
                   double *Xcosines, double *Ycosines, double *Zcosines) {
  int fin, varid, minc_icv, id, ndims, i, dim[256], length,
    minfound = 0, maxfound = 0;
  long start[3] = {0, 0, 0}, count[3];
  unsigned int plane1, plane2, line, shift[3], ct[3];
  nc_type type;
  char sign_type[MI_MAX_ATTSTR_LEN], name[MAX_NC_NAME];
  void *data;
  double range[2], swap;

  /* open minc file */
  fin = miopen((char *) filename, NC_NOWRITE);
  if(fin == MI_ERROR) return 0;

  /* get number of dimensions and type */
  id = ncvarid(fin, MIimage);
  (void) ncvarinq(fin, id, nullptr, &type, &ndims, dim, nullptr);
  if(ndims != 3) {
    fprintf(stderr, "unsupported %i dimensional minc file\n", ndims);
    return 0;
  }

  /* get sign */
  if ((miattgetstr(fin, id, MIsigntype, MI_MAX_ATTSTR_LEN, sign_type)
       == nullptr) || ((strcmp(sign_type, MI_UNSIGNED)!=0) &&
                    (strcmp(sign_type, MI_SIGNED)!=0))) {
    if (type == NC_BYTE)
      (void) strcpy(sign_type, MI_UNSIGNED);
    else
      (void) strcpy(sign_type, MI_SIGNED);
  }
  im->sign = (!strcmp(sign_type, MI_UNSIGNED)) ? SGN_UNSIGNED : SGN_SIGNED;

  switch(type) {
  case NC_CHAR:
  case NC_BYTE:
    im->wdim = 1;
    im->wordKind = WK_FIXED;
    im->sign = SGN_UNSIGNED;
    break;
  case NC_SHORT:
    im->wdim = 2;
    im->wordKind = WK_FIXED;
    break;
  case NC_LONG:
    im->wdim = 4;
    im->wordKind = WK_FIXED;
    break;
  case NC_FLOAT:
    im->wdim = 4;
    im->wordKind = WK_FLOAT;
    im->sign = SGN_UNKNOWN;
    break;
  case NC_DOUBLE:
    im->wdim = 8;
    im->wordKind = WK_FLOAT;
    im->sign = SGN_UNKNOWN;
    break;
  default:
    fprintf(stderr, "unsupported minc type\n");
    return 0;
    break;
  }


  /* get dimensions length and voxel size */
  ncopts = 0;
  for(i = 0; i < ndims; i++) {
    (void) ncdiminq(fin, dim[i], name, &count[i]);
    varid = ncvarid(fin, name);
    if(!strcmp(name, MIxspace)) {
      im->xdim = count[i];
      if(miattget1(fin, varid, MIstep, NC_DOUBLE, stepx) == MI_ERROR)
        *stepx = 1.0;
      im->vx = fabs(*stepx);
      if(miattget1(fin, varid, MIstart, NC_DOUBLE, startx) == MI_ERROR)
        *startx = 0.0;
      if(miattget(fin, varid, MIdirection_cosines, NC_DOUBLE, 3,
                  (void *) Xcosines, (int *) 0) == MI_ERROR ) {
        Xcosines[0] = 0.0;
        Xcosines[1] = 0.0;
        Xcosines[2] = 0.0;
      }
    }
    else if(!strcmp(name, MIyspace)) {
      im->ydim = count[i];
      if(miattget1(fin, varid, MIstep, NC_DOUBLE, stepy) == MI_ERROR)
        *stepy = 1.0;
      im->vy = fabs(*stepy);
      if(miattget1(fin, varid, MIstart, NC_DOUBLE, starty) == MI_ERROR)
        *starty = 0.0;
      if(miattget(fin, varid, MIdirection_cosines, NC_DOUBLE, 3,
                  (void *) Ycosines, (int *) 0) == MI_ERROR ) {
        Ycosines[0] = 0.0;
        Ycosines[1] = 0.0;
        Ycosines[2] = 0.0;
      }
    }
    else if(!strcmp(name, MIzspace)) {
      im->zdim = count[i];
      if(miattget1(fin, varid, MIstep, NC_DOUBLE, stepz) == MI_ERROR)
        *stepz = 1.0;
      im->vz = fabs(*stepz);
      if(miattget1(fin, varid, MIstart, NC_DOUBLE, startz) == MI_ERROR)
        *startz = 0.0;
      if(miattget(fin, varid, MIdirection_cosines, NC_DOUBLE, 3,
                  (void *) Zcosines, (int *) 0) == MI_ERROR ) {
        Zcosines[0] = 0.0;
        Zcosines[1] = 0.0;
        Zcosines[2] = 0.0;
      }
    } else {
      fprintf(stderr, "unsupported dimension %s\n", name);
      return 0;
    }
  }


  if( miattget( fin, id, MIvalid_range, NC_DOUBLE, 2,
                (void *) range, &length) == MI_ERROR || length != 2 ) {
    if( miattget1( fin, id, MIvalid_min, NC_DOUBLE,
                   (void *) &range[0] ) != MI_ERROR )
      minfound = 1;
    if( miattget1( fin, id, MIvalid_max, NC_DOUBLE,
                   (void *) &range[1] ) != MI_ERROR )
      maxfound = 1;
  }
  else {
    if( range[0] > range[1] ) {
      swap = range[0];
      range[0] = range[1];
      range[1] = swap;
    }
    minfound = 1;
    maxfound = 1;
  }

  if(!minfound) {
    switch(type) {
    case NC_CHAR:
    case NC_BYTE:
      range[0] = 0;
      break;
    case NC_SHORT:
      if(im->sign == SGN_SIGNED) range[0] = -32768.0;
      else range[0] = 0;
      break;
    case NC_LONG:
      if(im->sign == SGN_SIGNED) range[0] = -2147483648.0;
      else range[0] = 0;
      break;
    default:
      range[0] = 0.0;
      break;
    }
  }
  if(!maxfound) {
    switch(type) {
    case NC_CHAR:
    case NC_BYTE:
      range[1] = 255.0;
      break;
    case NC_SHORT:
      if(im->sign == SGN_SIGNED) range[1] = 32767.0;
      else range[1] = 65535.0;
      break;
    case NC_LONG:
      if(im->sign == SGN_SIGNED) range[1] = 2147483647.0;
      else range[1] = 4294967295.0;
      break;
    default:
      range[1] = 1.0;
      break;
    }
  }

  ncopts = NC_VERBOSE | NC_FATAL;
  im->vdim = 1;

  /* read data bytes */
  im->data = ImageIO_alloc(im->xdim * im->ydim * im->zdim * im->vdim *
                           im->wdim);

  minc_icv = miicv_create();
  (void) miicv_setint( minc_icv, MI_ICV_TYPE, type );
  (void) miicv_setstr( minc_icv, MI_ICV_SIGN, sign_type );
  (void) miicv_setint( minc_icv, MI_ICV_DO_NORM, 1 );
  (void) miicv_setint( minc_icv, MI_ICV_DO_FILLVALUE, 1 );
  (void) miicv_setdbl( minc_icv, MI_ICV_VALID_MIN, range[0] );
  (void) miicv_setdbl( minc_icv, MI_ICV_VALID_MAX, range[1] );

  (void) miicv_attach( minc_icv, fin, id );

  if(miicv_get(minc_icv, start, count, im->data) ==  MI_ERROR) {
    fprintf(stderr, "error while reading file\n");
    return 0;
  }
  (void) miicv_detach( minc_icv );
  (void) miicv_free( minc_icv );

  ImageIO_closeImage(im);

  /* order data in ZYX */
  (void) ncdiminq(fin, dim[0], name, nullptr);
  if(!strcmp(name, MIzspace)) {
    (void) ncdiminq(fin, dim[1], name, nullptr);
    /* file is ZYX */
    if(!strcmp(name, MIyspace)) {
      miclose(fin);
      return 1;
    }
  }

  (void) ncdiminq(fin, dim[0], name, nullptr);
  /* file is ZXY */
  if(!strcmp(name, MIzspace)) {
    shift[0] = 0;
    shift[1] = 2;
    shift[2] = 1;
    plane2 = im->xdim * im->ydim;
    line = im->ydim;
  }
  /* file is Y first */
  else if(!strcmp(name, MIyspace)) {
    shift[0] = 1;
    plane2 = im->xdim * im->zdim;
    (void) ncdiminq(fin, dim[1], name, nullptr);
    /* file is YXZ */
    if(!strcmp(name, MIxspace)) {
      shift[1] = 2;
      shift[2] = 0;
      line = im->zdim;
    }
    /* file is YZX */
    else {
      shift[1] = 0;
      shift[2] = 2;
      line = im->xdim;
    }
  }
  /* file is X first */
  else {
    shift[0] = 2;
    plane2 = im->ydim * im->zdim;
    (void) ncdiminq(fin, dim[1], name, nullptr);
    /* file is XYZ */
    if(!strcmp(name, MIyspace)) {
      shift[1] = 1;
      shift[2] = 0;
      line = im->zdim;
    }
    /* file is XZY */
    else {
      shift[1] = 0;
      shift[2] = 1;
      line = im->ydim;
    }
  }

  data = ImageIO_alloc(im->xdim * im->ydim * im->zdim * im->vdim * im->wdim);
  plane1 = im->xdim * im->ydim;

  for(ct[0] = 0; ct[0] < im->zdim; ct[0]++) {
    for(ct[1] = 0; ct[1] < im->ydim; ct[1]++) {
      for(ct[2] = 0; ct[2] < im->xdim; ct[2]++) {
        memcpy((unsigned char *) data +
               (ct[0] * plane1 + ct[1] * im->xdim + ct[2]) * im->wdim,
               (unsigned char *) im->data +
               (ct[shift[0]] * plane2 + ct[shift[1]] * line + ct[shift[2]])
               * im->wdim, im->wdim);
      }
    }
  }

  ImageIO_free(im->data);
  im->data = data;

  /* close file */
  miclose(fin);
  return 1;
}



CGAL_INLINE_FUNCTION
int writeMincFile( const _image* im, const char *filename,
                   const char *sourcename,
                   double startx, double starty, double startz,
                   double stepx, double stepy, double stepz,
                   const double *Xcosine, const double *Ycosine,
                   const double *Zcosine, const double *range ) {
  int i, cdfid, minc_icv, dim_vars[3], dim_ids[3], img_var_id, min_id, max_id,
    src_cdfid, src_img_var, varid, n_excluded, excluded_vars[10],
    src_min_id, src_max_id, src_root_id, old_att_length, same_file = 0;
  long org[3] = {0, 0, 0}, count[3] = { im->zdim, im->ydim, im->xdim };
  char *dim_names[3] = { MIzspace, MIyspace, MIxspace }, *history,
       *excluded_list[3] = { MIxspace, MIyspace, MIzspace }, *newname;
  nc_type type, datatype;
  double start[3] = { startz, starty, startx },
            vx[3] = { stepz, stepy, stepx };

  /* if filename is the same file that sourcename, it is needed to create
     a temporary file for output */
  newname = (char *) filename;
  if(sourcename) {
    struct stat out, org;
    same_file = (stat(sourcename, &org) == 0) && (stat(filename, &out) == 0)
      && (org.st_ino == out.st_ino);
    if(same_file) {
      newname = (char *) malloc((strlen(filename) + 10) * sizeof(char));
      for(i = strlen(filename) - 1; i > 0 && filename[i] != '/'; i--) ;
      if(filename[i] == '/') {
        strncpy(newname, filename, i + 1);
        newname[i+1] = '\0';
        strcat(newname, "#TMP#");
        strcat(newname, filename + i + 1);
      }
      else
        sprintf(newname, "#TMP#%s", filename);
    }
  }

  /* open file */
  cdfid = micreate( (char *) newname, NC_CLOBBER );
  if( cdfid == MI_ERROR ) {
    fprintf(stderr, "Error: opening MINC file \"%s\".\n", filename );
    return -1;
  }

  /* set word type */
  switch(im->wordKind) {
  case WK_FIXED:
    if(im->wdim == 1) type = NC_BYTE;
    else if(im->wdim == 2) type = NC_SHORT;
    else type = NC_LONG;
    break;
  case WK_FLOAT:
  default:
    if(im->wdim == 4) type = NC_FLOAT;
    else type = NC_DOUBLE;
    break;
  }

  /* create dimensions variables */
  for(i = 0; i < 3; i++)
    dim_vars[i] = ncdimdef( cdfid, dim_names[i], count[i] );

  for(i = 0; i < 3; i++ ) {
    dim_ids[i] = micreate_std_variable( cdfid, dim_names[i], NC_DOUBLE,
                                        0, nullptr);
    if( dim_ids[i] < 0 ) return -1;
    (void) miattputdbl( cdfid, dim_ids[i], MIstep, vx[i]);
    (void) miattputdbl( cdfid, dim_ids[i], MIstart, start[i]);
  }
  if(Zcosine[0] != 0.0 || Zcosine[1] != 0.0 || Zcosine[2] != 0.0)
    (void) ncattput(cdfid, dim_ids[0], MIdirection_cosines, NC_DOUBLE, 3,
                    Zcosine);
  if(Ycosine[0] != 0.0 || Ycosine[1] != 0.0 || Ycosine[2] != 0.0)
    (void) ncattput(cdfid, dim_ids[1], MIdirection_cosines, NC_DOUBLE, 3,
                    Ycosine);
  if(Xcosine[0] != 0.0 || Xcosine[1] != 0.0 || Xcosine[2] != 0.0)
    (void) ncattput(cdfid, dim_ids[2], MIdirection_cosines, NC_DOUBLE, 3,
                    Xcosine);

  /* create image variable */
  img_var_id = micreate_std_variable( cdfid, MIimage, type, 3, dim_vars );
  (void) miattputstr( cdfid, img_var_id, MIcomplete, MI_FALSE );
  if(im->sign == SGN_SIGNED)
    (void) miattputstr( cdfid, img_var_id, MIsigntype, MI_SIGNED );
  else
    (void) miattputstr( cdfid, img_var_id, MIsigntype, MI_UNSIGNED );
  (void) ncattput(cdfid, img_var_id, MIvalid_range ,NC_DOUBLE, 2, range);

  min_id = micreate_std_variable(cdfid, MIimagemin, NC_DOUBLE, 0, 0 );
  max_id = micreate_std_variable(cdfid, MIimagemax, NC_DOUBLE, 0, 0 );

  /* if the image comes from a minc file, duplicate all minc attributes */
  if(sourcename) {
    if( (src_cdfid = miopen((char *) sourcename, NC_NOWRITE)) == MI_ERROR ) {
      fprintf(stderr, "Error opening %s\n", sourcename );
      return -1;
    }

    n_excluded = 0;
    for( i = 0; i < 3; i++ ) {
      if( (varid = ncvarid(src_cdfid, excluded_list[i] )) != MI_ERROR )
        excluded_vars[n_excluded++] = varid;
    }

    if( (src_img_var = ncvarid(src_cdfid, MIimage )) != MI_ERROR )
      excluded_vars[n_excluded++] = src_img_var;
    if( (src_max_id = ncvarid(src_cdfid, MIimagemax )) != MI_ERROR )
      excluded_vars[n_excluded++] = src_max_id;
    if( (src_min_id = ncvarid(src_cdfid, MIimagemin )) != MI_ERROR )
      excluded_vars[n_excluded++] = src_min_id;
    if( (src_root_id = ncvarid(src_cdfid, MIrootvariable )) != MI_ERROR )
      excluded_vars[n_excluded++] = src_root_id;

    (void) micopy_all_var_defs( src_cdfid, cdfid, n_excluded, excluded_vars );

    if( src_img_var != MI_ERROR )
      (void) micopy_all_atts( src_cdfid, src_img_var, cdfid, img_var_id );
    if( src_root_id != MI_ERROR )
      (void) micopy_all_atts( src_cdfid, src_root_id, cdfid,
                              ncvarid( cdfid, MIrootvariable) );

    ncopts = 0;
    if(ncattinq(cdfid, NC_GLOBAL, MIhistory, &datatype, &old_att_length)
       == MI_ERROR || datatype != NC_CHAR )
      old_att_length = 0;

    history = (char *) malloc(sizeof(char) * old_att_length);
    history[0] = (char) 0;
    (void) miattgetstr( cdfid, NC_GLOBAL, MIhistory, old_att_length+1,
                        history );
    ncopts = NC_VERBOSE | NC_FATAL;
    (void) miattputstr( cdfid, NC_GLOBAL, MIhistory, history );
    free(history);

    (void) miclose( src_cdfid );
  }


  /* write data */
  if(ncendef(cdfid) == MI_ERROR) return -1;
  minc_icv = miicv_create();
  (void) miicv_setint( minc_icv, MI_ICV_TYPE, type);
  (void) miicv_setstr( minc_icv, MI_ICV_SIGN,
                       (im->sign == SGN_SIGNED) ? MI_SIGNED : MI_UNSIGNED );
  (void) miicv_setint( minc_icv, MI_ICV_DO_NORM, 1 );
  (void) miicv_setint( minc_icv, MI_ICV_USER_NORM, 1 );
  (void) miicv_attach( minc_icv, cdfid, img_var_id );

  if( miicv_put(minc_icv, org, count, im->data) == MI_ERROR )
    return -1;

  (void) miattputstr( cdfid, img_var_id, MIcomplete, MI_TRUE);
  (void) miclose( cdfid );
  (void) miicv_free( minc_icv );

  if(same_file) {
    unlink(filename);
    link(newname, filename);
    unlink(newname);
    free(newname);
  }

  return 0;
}


#endif
