	#.file	1 "workaround_4_irix.c"
	.option pic2
	.section	.text
	.text
	.align	2
	.globl	CGAL_workaround_IRIX_set_FPU_cw
	.ent	CGAL_workaround_IRIX_set_FPU_cw
CGAL_workaround_IRIX_set_FPU_cw:
	.frame	$fp,48,$31		# vars= 16, regs= 2/0, args= 0, extra= 16
	.mask	0x50000000,-8
	.fmask	0x00000000,0
	subu	$sp,$sp,48
	sd	$fp,40($sp)
	sd	$28,32($sp)
	move	$fp,$sp
	.set	noat
	lui	$1,%hi(%neg(%gp_rel(CGAL_workaround_IRIX_set_FPU_cw)))
	addiu	$1,$1,%lo(%neg(%gp_rel(CGAL_workaround_IRIX_set_FPU_cw)))
	daddu	$gp,$1,$25
	.set	at
	sw	$4,16($fp)
	lw	$2,16($fp)
 #APP
	ctc1 $2,$31
 #NO_APP
.L1:
	move	$sp,$fp
	ld	$fp,40($sp)
	ld	$28,32($sp)
	addu	$sp,$sp,48
	j	$31
	.end	CGAL_workaround_IRIX_set_FPU_cw
	.align	2
	.globl	CGAL_workaround_IRIX_get_FPU_cw
	.ent	CGAL_workaround_IRIX_get_FPU_cw
CGAL_workaround_IRIX_get_FPU_cw:
	.frame	$fp,48,$31		# vars= 16, regs= 2/0, args= 0, extra= 16
	.mask	0x50000000,-8
	.fmask	0x00000000,0
	subu	$sp,$sp,48
	sd	$fp,40($sp)
	sd	$28,32($sp)
	move	$fp,$sp
	.set	noat
	lui	$1,%hi(%neg(%gp_rel(CGAL_workaround_IRIX_get_FPU_cw)))
	addiu	$1,$1,%lo(%neg(%gp_rel(CGAL_workaround_IRIX_get_FPU_cw)))
	daddu	$gp,$1,$25
	.set	at
 #APP
	cfc1 $4,$31
 #NO_APP
	sw	$4,16($fp)
	lw	$3,16($fp)
	move	$2,$3
	b	.L2
.L2:
	move	$sp,$fp
	ld	$fp,40($sp)
	ld	$28,32($sp)
	addu	$sp,$sp,48
	j	$31
	.end	CGAL_workaround_IRIX_get_FPU_cw
