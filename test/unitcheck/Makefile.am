
if HAVE_CHECKMK

bin_PROGRAMS=lctest lcutiltest

lctest_SOURCES=basic_complete.c
lctest_CPPFLAGS=-I../../libqs/include -Wall @CHECK_CFLAGS@ -I../../libqsutil -g
lctest_LDADD =../../libqs/libqsearch.la ../../libqsutil/libqsutil.la -lm @CHECK_LIBS@

lcutiltest_SOURCES=basic_qsutil.c
lcutiltest_CPPFLAGS=-I../../libqs/include -Wall @CHECK_CFLAGS@ -I../../libqsutil -g
lcutiltest_LDADD =../../libqs/libqsearch.la ../../libqsutil/libqsutil.la -lm @CHECK_LIBS@

.ts.c:
	checkmk $< >$@

endif

EXTRA_DIST=basic_complete.ts \
  basic_qsutil.ts

clean-local:
	rm -f *.c
