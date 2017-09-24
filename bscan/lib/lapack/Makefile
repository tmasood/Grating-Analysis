#######################################################################
#  This is the makefile to create a library for LAPACK.
#  The files are organized as follows:
#       ALLAUX -- Auxiliary routines called from all precisions
#       SCLAUX -- Auxiliary routines called from both REAL and COMPLEX
#       DZLAUX -- Auxiliary routines called from both DOUBLE PRECISION
#                 and COMPLEX*16
#       SLASRC -- Single precision real LAPACK routines
#       CLASRC -- Single precision complex LAPACK routines
#       DLASRC -- Double precision real LAPACK routines
#       ZLASRC -- Double precision complex LAPACK routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       lapack.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#  ***Note***
#  The functions lsame, second, dsecnd, slamch, and dlamch may have
#  to be installed before compiling the library.  Refer to the
#  installation guide, LAPACK Working Note 41, for instructions.
#
#----------------------------------------------------------------------
#
#  Edward Anderson, University of Tennessee
#  March 26, 1990. 
#  Susan Ostrouchov, Last updated September 30, 1994
#
#######################################################################
CC		= gcc
CFLAGS		= -fPIC -c -O2
#######################################################################

ALLAUX = ilaenv.o lsame.o lsamen.o xerbla.o

DZLAUX = \
   dbdsqr.o ddisna.o dlabad.o dlacpy.o dladiv.o dlae2.o  dlaebz.o \
   dlaed0.o dlaed1.o dlaed2.o dlaed3.o dlaed4.o dlaed5.o dlaed6.o \
   dlaed7.o dlaed8.o dlaed9.o dlaeda.o dlaev2.o dlagtf.o \
   dlagts.o dlamch.o dlamrg.o dlanst.o dlapy2.o dlapy3.o dlarnv.o \
   dlartg.o dlaruv.o \
   dlas2.o  dlascl.o dlaset.o dlasq1.o dlasq2.o dlasq3.o dlasq4.o \
   dlasr.o  dlasrt.o dlassq.o dlasv2.o dpttrf.o dstebz.o dstedc.o \
   dsteqr.o dsterf.o dsecnd.o

DLASRC = \
   dgbbrd.o dgbcon.o dgbequ.o dgbrfs.o dgbsv.o  \
   dgbsvx.o dgbtf2.o dgbtrf.o dgbtrs.o dgebak.o dgebal.o dgebd2.o \
   dgebrd.o dgecon.o dgeequ.o dgees.o  dgeesx.o dgeev.o  dgeevx.o \
   dgegs.o  dgegv.o  dgehd2.o dgehrd.o dgelq2.o dgelqf.o \
   dgels.o  dgelss.o dgelsx.o dgeql2.o dgeqlf.o \
   dgeqpf.o dgeqr2.o dgeqrf.o dgerfs.o dgerq2.o dgerqf.o \
   dgesv.o  dgesvd.o dgesvx.o dgetf2.o dgetrf.o dgetri.o \
   dgetrs.o dggbak.o dggbal.o dggglm.o dgghrd.o dgglse.o dggqrf.o \
   dggrqf.o dggsvd.o dggsvp.o dgtcon.o dgtrfs.o dgtsv.o  \
   dgtsvx.o dgttrf.o dgttrs.o dhgeqz.o \
   dhsein.o dhseqr.o dlabrd.o dlacon.o \
   dlaein.o dlaexc.o dlag2.o  dlags2.o dlagtm.o dlahqr.o \
   dlahrd.o dlaic1.o dlaln2.o \
   dlangb.o dlange.o dlangt.o dlanhs.o dlansb.o dlansp.o \
   dlansy.o dlantb.o dlantp.o dlantr.o dlanv2.o \
   dlapll.o dlapmt.o \
   dlaqgb.o dlaqge.o dlaqsb.o dlaqsp.o dlaqsy.o dlaqtr.o dlar2v.o \
   dlarf.o  dlarfb.o dlarfg.o dlarft.o dlarfx.o dlargv.o dlartv.o \
   dlaswp.o dlasy2.o dlasyf.o \
   dlatbs.o dlatps.o dlatrd.o dlatrs.o dlatzm.o dlauu2.o dlauum.o \
   dopgtr.o dopmtr.o dorg2l.o dorg2r.o \
   dorgbr.o dorghr.o dorgl2.o dorglq.o dorgql.o dorgqr.o dorgr2.o \
   dorgrq.o dorgtr.o dorm2l.o dorm2r.o \
   dormbr.o dormhr.o dorml2.o dormlq.o dormql.o dormqr.o dormr2.o \
   dormrq.o dormtr.o dpbcon.o dpbequ.o dpbrfs.o dpbstf.o \
   dpbsv.o  dpbsvx.o \
   dpbtf2.o dpbtrf.o dpbtrs.o dpocon.o dpoequ.o dporfs.o dposv.o  \
   dposvx.o dpotf2.o dpotrf.o dpotri.o dpotrs.o dppcon.o dppequ.o \
   dpprfs.o dppsv.o  dppsvx.o dpptrf.o dpptri.o dpptrs.o dptcon.o \
   dpteqr.o dptrfs.o dptsv.o  dptsvx.o dpttrs.o drscl.o  \
   dsbev.o  dsbevd.o dsbevx.o dsbgst.o dsbgv.o  dsbtrd.o dspcon.o \
   dspev.o  dspevd.o dspevx.o dspgst.o \
   dspgv.o  dsprfs.o dspsv.o  dspsvx.o dsptrd.o dsptrf.o dsptri.o \
   dsptrs.o dstein.o dstev.o  dstevd.o dstevx.o \
   dsycon.o dsyev.o  dsyevd.o \
   dsyevx.o dsygs2.o dsygst.o dsygv.o  dsyrfs.o dsysv.o  dsysvx.o \
   dsytd2.o dsytf2.o dsytrd.o dsytrf.o dsytri.o dsytrs.o dtbcon.o \
   dtbrfs.o dtbtrs.o dtgevc.o dtgsja.o dtpcon.o dtprfs.o dtptri.o \
   dtptrs.o \
   dtrcon.o dtrevc.o dtrexc.o dtrrfs.o dtrsen.o dtrsna.o dtrsyl.o \
   dtrti2.o dtrtri.o dtrtrs.o dtzrqf.o

ZLASRC = \
   zbdsqr.o zgbbrd.o zgbcon.o zgbequ.o zgbrfs.o zgbsv.o  zgbsvx.o \
   zgbtf2.o zgbtrf.o zgbtrs.o zgebak.o zgebal.o zgebd2.o zgebrd.o \
   zgecon.o zgeequ.o zgees.o  zgeesx.o zgeev.o  zgeevx.o \
   zgegs.o  zgegv.o  zgehd2.o zgehrd.o zgelq2.o zgelqf.o \
   zgels.o  zgelss.o zgelsx.o zgeql2.o zgeqlf.o zgeqpf.o \
   zgeqr2.o zgeqrf.o zgerfs.o zgerq2.o zgerqf.o \
   zgesv.o  zgesvd.o zgesvx.o zgetf2.o zgetrf.o zgetri.o zgetrs.o \
   zggbak.o zggbal.o zggglm.o zgghrd.o zgglse.o zggqrf.o zggrqf.o \
   zggsvd.o zggsvp.o \
   zgtcon.o zgtrfs.o zgtsv.o  zgtsvx.o zgttrf.o zgttrs.o zhbev.o  \
   zhbevd.o zhbevx.o zhbgst.o zhbgv.o  zhbtrd.o zhecon.o zheev.o  \
   zheevd.o zheevx.o zhegs2.o zhegst.o \
   zhegv.o  zherfs.o zhesv.o  zhesvx.o zhetd2.o zhetf2.o zhetrd.o \
   zhetrf.o zhetri.o zhetrs.o zhgeqz.o zhpcon.o zhpev.o  zhpevd.o \
   zhpevx.o zhpgst.o zhpgv.o  zhprfs.o zhpsv.o  zhpsvx.o \
   zhptrd.o zhptrf.o zhptri.o zhptrs.o zhsein.o zhseqr.o zlabrd.o \
   zlacgv.o zlacon.o zlacpy.o zlacrm.o zlacrt.o zladiv.o \
   zlaed0.o zlaed7.o zlaed8.o \
   zlaein.o zlaesy.o zlaev2.o zlags2.o zlagtm.o \
   zlahef.o zlahqr.o \
   zlahrd.o zlaic1.o zlangb.o zlange.o zlangt.o zlanhb.o zlanhe.o \
   zlanhp.o zlanhs.o zlanht.o zlansb.o zlansp.o zlansy.o zlantb.o \
   zlantp.o zlantr.o zlapll.o zlapmt.o zlaqgb.o zlaqge.o \
   zlaqhb.o zlaqhe.o zlaqhp.o zlaqsb.o \
   zlaqsp.o zlaqsy.o zlar2v.o zlarf.o  zlarfb.o zlarfg.o zlarft.o \
   zlarfx.o zlargv.o zlarnv.o zlartg.o zlartv.o \
   zlascl.o zlaset.o zlasr.o  zlassq.o zlaswp.o zlasyf.o \
   zlatbs.o zlatps.o zlatrd.o zlatrs.o zlatzm.o zlauu2.o zlauum.o \
   zpbcon.o zpbequ.o zpbrfs.o zpbstf.o zpbsv.o  \
   zpbsvx.o zpbtf2.o zpbtrf.o zpbtrs.o zpocon.o zpoequ.o zporfs.o \
   zposv.o  zposvx.o zpotf2.o zpotrf.o zpotri.o zpotrs.o zppcon.o \
   zppequ.o zpprfs.o zppsv.o  zppsvx.o zpptrf.o zpptri.o zpptrs.o \
   zptcon.o zpteqr.o zptrfs.o zptsv.o  zptsvx.o zpttrf.o zpttrs.o \
   zrot.o   zspcon.o zspmv.o  zspr.o   zsprfs.o zspsv.o  \
   zspsvx.o zsptrf.o zsptri.o zsptrs.o zdrot.o  zdrscl.o zstedc.o \
   zstein.o zsteqr.o zsycon.o zsymv.o  \
   zsyr.o   zsyrfs.o zsysv.o  zsysvx.o zsytf2.o zsytrf.o zsytri.o \
   zsytrs.o ztbcon.o ztbrfs.o ztbtrs.o ztgevc.o ztgsja.o ztpcon.o \
   ztprfs.o ztptri.o \
   ztptrs.o ztrcon.o ztrevc.o ztrexc.o ztrrfs.o ztrsen.o ztrsna.o \
   ztrsyl.o ztrti2.o ztrtri.o ztrtrs.o ztzrqf.o zung2l.o zung2r.o \
   zungbr.o zunghr.o zungl2.o zunglq.o zungql.o zungqr.o zungr2.o \
   zungrq.o zungtr.o zunm2l.o zunm2r.o zunmbr.o zunmhr.o zunml2.o \
   zunmlq.o zunmql.o zunmqr.o zunmr2.o zunmrq.o zunmtr.o zupgtr.o \
   zupmtr.o izmax1.o dzsum1.o

library	: $(ALLAUX) $(DZLAUX) $(DLASRC) $(ZLASRC)
	gcc -shared -Wl,-soname,liblapack.so.1 -o liblapack.so.1.0\
	$(ALLAUX) $(DZLAUX) $(DLASRC) $(ZLASRC)

install : 
	mv liblapack.so.1.0 ../
	ln -sf ../liblapack.so.1.0 ../liblapack.so.1
	ln -sf ../liblapack.so.1 ../liblapack.so

clean:
	rm -f *.o

.c.o: 
	$(CC) $(CFLAGS) -c $*.c















