if HAVE_LD_VERSION_SCRIPT
    VERSION_LDFLAGS = -Wl,--version-script=$(top_srcdir)/lib/qsutil.map
else
    VERSION_LDFLAGS = -export-symbols $(top_srcdir)/lib/qsutil.sym
endif

lib_LTLIBRARIES = libqsutil.la
libqsutil_la_SOURCES = qsutil.c qsutil.h

libqsutil_la_CPPFLAGS = -I$(top_srcdir)/libqs/include -Wall
libqsutil_la_LDFLAGS = $(VERSION_LDFLAGS)
