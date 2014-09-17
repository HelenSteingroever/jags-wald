noinst_LTLIBRARIES = walddist.la

walddist_la_CPPFLAGS = -I$(top_srcdir)/src \
												 -I$(includedir)/JAGS

walddist_la_LDFLAGS = -no-undefined -module -avoid-version

walddist_la_SOURCES = DWald.cc DWaldGamma.cc

noinst_HEADERS = DWald.h DWaldGamma.h
