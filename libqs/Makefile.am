
if HAVE_LD_VERSION_SCRIPT
    VERSION_LDFLAGS = -Wl,--version-script=$(top_srcdir)/lib/qsearch.map
else
    VERSION_LDFLAGS = -export-symbols $(top_srcdir)/lib/qsearch.sym
endif

lib_LTLIBRARIES = libqsearch.la
libqsearch_la_SOURCES = quartet_tree.c libqs.c mcmc.c
libqsearch_la_CPPFLAGS = -I$(top_srcdir)/include -Wall -O3
libqsearch_la_CFLAGS = -I$(top_srcdir)/include -Wall -O3
libqsearch_la_LDFLAGS = $(VERSION_LDFLAGS) -O3
libqsearch_la_LIBADD =
include_HEADERS =   include/qsearch.h

pubincludedir = $(includedir)/qsearch
pubinclude_HEADERS =                                              \
                       include/qsearch/quartet_tree.h             \
                       include/qsearch/libqs.h
