ACLOCAL_AMFLAGS = -I m4
SUBDIRS = libqs libqsutil tools test

EXTRA_DIST = lib/qsearch.map lib/qsearch.sym lib/qsutil.map lib/qsutil.sym \
             lib/libqsearch.pc.in
dist_man_MANS = man/maketree.1

pkgconfigdir = $(libdir)/pkgconfig
nodist_pkgconfig_DATA = lib/libqsearch.pc

lib/libqsearch.pc: $(srcdir)/lib/libqsearch.pc.in
	sed -e 's![@]prefix[@]!$(prefix)!g' \
			-e 's![@]exec_prefix[@]!$(exec_prefix)!g' \
			-e 's![@]includedir[@]!$(includedir)!g' \
			-e 's![@]libdir[@]!$(libdir)!g' \
			-e 's![@]PACKAGE_VERSION[@]!$(PACKAGE_VERSION)!g' \
			-e 's![@]QSEARCH1_LDFLAGS[@]!$(QSEARCH1_LDFLAGS)!g' \
			-e 's![@]QSEARCH1_LIBS[@]!$(QSEARCH1_LIBS)!g' \
			$(srcdir)/lib/libqsearch.pc.in > $@
