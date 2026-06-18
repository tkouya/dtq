LDADD = $(top_builddir)/src/libqd.la
AM_CPPFLAGS = -I$(top_builddir) -I$(top_builddir)/include -I$(top_srcdir)/include

TESTS = qd_test pslq_test c_test
check_PROGRAMS = qd_test pslq_test c_test
EXTRA_PROGRAMS = qd_timer quadt_test huge

dist_noinst_DATA = coeff.dat

CLEANFILES=qd_timer quadt_test huge

qd_test_SOURCES = qd_test.cpp
pslq_test_SOURCES = pslq.h pslq_test.cpp tictoc.cpp tictoc.h
qd_timer_SOURCES = qd_timer.cpp tictoc.cpp tictoc.h
quadt_test_SOURCES = quadt_test.cpp tictoc.cpp quadt.h tictoc.h
huge_SOURCES = huge.cpp

c_test_SOURCES = c_test.c
c_test_LINK = $(CXXLINK)

time: qd_timer$(EXEEXT)
	./qd_timer$(EXEEXT)

demo: $(check_PROGRAMS) $(EXTRA_PROGRAMS)

.PHONY: time demo
