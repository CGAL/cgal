	#.file	1 "workaround_4_irix.c"
	.set	nobopt
	.option pic2
	.section .text,0x1,0x6,4,8
	.text
	.align	2
	.align	3
	.globl	CGAL_workaround_IRIX_set_FPU_cw
	.ent	CGAL_workaround_IRIX_set_FPU_cw
CGAL_workaround_IRIX_set_FPU_cw:
.LFB1:
	.frame	$sp,32,$31		# vars= 0, regs= 1/0, args= 0, extra= 16
	.mask	0x10000000,-16
	.fmask	0x00000000,0
	dsubu	$sp,$sp,32
.LCFI0:
	sd	$28,16($sp)
.LCFI1:
	.set	noat
	lui	$1,%hi(%neg(%gp_rel(CGAL_workaround_IRIX_set_FPU_cw)))
	addiu	$1,$1,%lo(%neg(%gp_rel(CGAL_workaround_IRIX_set_FPU_cw)))
	daddu	$gp,$1,$25
	.set	at
 #APP
	ctc1 $4,$31
 #NO_APP
	ld	$28,16($sp)
	#nop
	.set	noreorder
	.set	nomacro
	j	$31
	daddu	$sp,$sp,32
	.set	macro
	.set	reorder

.LFE1:
	.end	CGAL_workaround_IRIX_set_FPU_cw
	.align	2
	.align	3
	.globl	CGAL_workaround_IRIX_get_FPU_cw
	.ent	CGAL_workaround_IRIX_get_FPU_cw
CGAL_workaround_IRIX_get_FPU_cw:
.LFB2:
	.frame	$sp,32,$31		# vars= 0, regs= 1/0, args= 0, extra= 16
	.mask	0x10000000,-16
	.fmask	0x00000000,0
	dsubu	$sp,$sp,32
.LCFI2:
	sd	$28,16($sp)
.LCFI3:
	.set	noat
	lui	$1,%hi(%neg(%gp_rel(CGAL_workaround_IRIX_get_FPU_cw)))
	addiu	$1,$1,%lo(%neg(%gp_rel(CGAL_workaround_IRIX_get_FPU_cw)))
	daddu	$gp,$1,$25
	.set	at
 #APP
	cfc1 $2,$31
 #NO_APP
	ld	$28,16($sp)
	#nop
	.set	noreorder
	.set	nomacro
	j	$31
	daddu	$sp,$sp,32
	.set	macro
	.set	reorder

.LFE2:
	.end	CGAL_workaround_IRIX_get_FPU_cw
