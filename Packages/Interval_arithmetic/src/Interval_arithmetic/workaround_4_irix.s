	.file	1 "Interval_arithmetic/workaround_4_irix.c"
	.option pic2
gcc2_compiled.:
__gnu_compiled_c:
	.text
	.align	2
	.globl	CGAL_workaround_IRIX_set_FPU_cw
	.ent	CGAL_workaround_IRIX_set_FPU_cw
CGAL_workaround_IRIX_set_FPU_cw:
	.frame	$fp,16,$31		# vars= 0, regs= 2/0, args= 0, extra= 8
	.mask	0x50000000,-4
	.fmask	0x00000000,0
	.set	noreorder
	.cpload	$25
	.set	reorder
	subu	$sp,$sp,16
	.cprestore 0
	sw	$fp,12($sp)
	sw	$28,8($sp)
	move	$fp,$sp
	sw	$4,16($fp)
	lw	$2,16($fp)
 #APP
	ctc1 $2,$31
 #NO_APP
$L1:
	move	$sp,$fp
	lw	$fp,12($sp)
	addu	$sp,$sp,16
	j	$31
	.end	CGAL_workaround_IRIX_set_FPU_cw
	.align	2
	.globl	CGAL_workaround_IRIX_get_FPU_cw
	.ent	CGAL_workaround_IRIX_get_FPU_cw
CGAL_workaround_IRIX_get_FPU_cw:
	.frame	$fp,24,$31		# vars= 8, regs= 2/0, args= 0, extra= 8
	.mask	0x50000000,-4
	.fmask	0x00000000,0
	.set	noreorder
	.cpload	$25
	.set	reorder
	subu	$sp,$sp,24
	.cprestore 0
	sw	$fp,20($sp)
	sw	$28,16($sp)
	move	$fp,$sp
 #APP
	cfc1 $4,$31
 #NO_APP
	sw	$4,8($fp)
	lw	$3,8($fp)
	move	$2,$3
	j	$L2
$L2:
	move	$sp,$fp
	lw	$fp,20($sp)
	addu	$sp,$sp,24
	j	$31
	.end	CGAL_workaround_IRIX_get_FPU_cw
