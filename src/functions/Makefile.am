noinst_LTLIBRARIES = waldfunc.la

waldfunc_la_CPPFLAGS = -I$(top_srcdir)/src \
                       -I$(includedir)/JAGS

waldfunc_la_LDFLAGS = -no-undefined -module -avoid-version

waldfunc_la_SOURCES = DSWaldfunc.cc PSWaldfunc.cc

noinst_HEADERS = DSWaldfunc.h PSWaldfunc.h
