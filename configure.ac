AC_PREREQ([2.69])
AC_INIT([libqsearch], [1.0.9], [rudi@morefuntechnologies.com])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_LIBOBJ_DIR([portable])
AC_CONFIG_MACRO_DIR([m4])
AM_INIT_AUTOMAKE()
LT_INIT()
LT_PREREQ([2.2])
AC_PROG_CC()
AC_CONFIG_SRCDIR([tools/maketree.c])
AC_CONFIG_HEADERS([config.h])
AC_PROG_CC
RRA_LD_VERSION_SCRIPT
AC_CONFIG_FILES([Makefile
                 libqs/Makefile
                 libqsutil/Makefile
                 tools/Makefile
                 test/Makefile
                 test/unitcheck/Makefile
                ])
AC_CHECK_HEADER([getopt.h], [ flag_header_getopt_found=1 ])

if test "x$flag_header_getopt_found" = "x" ; then
  AC_ERROR([Must have getopt.h header])
fi

AC_PATH_PROG([HAVE_CHECKMK_PATH], [checkmk], [notfound])

AM_CONDITIONAL([HAVE_CHECKMK], [ test ! "x$HAVE_CHECKMK_PATH" = "xnotfound" ])

if test "x$HAVE_CHECKMK_PATH" = "xnotfound" ; then
  AC_WARN([checkmk not found please install "check" package for unit tests])
else
  CHECK_CFLAGS=$(pkg-config check --cflags)
  CHECK_LIBS=$(pkg-config check --libs)
fi

AC_SUBST([CHECK_CFLAGS])
AC_SUBST([CHECK_LIBS])

AC_OUTPUT
