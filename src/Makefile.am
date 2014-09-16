SUBDIRS = distributions

jagsmod_LTLIBRARIES = wald.la

wald_la_SOURCES = wald.cc

wald_la_CPPFLAGS = -I$(includedir)/JAGS

wald_la_LIBADD = distributions/walddist.la

if WINDOWS
wald_la_LIBADD += -ljags-$(JAGS_MAJOR) -ljrmath-0
else
wald_la_LIBADD += -ljags -ljrmath
endif

wald_la_LDFLAGS = -module -avoid-version
if WINDOWS
wald_la_LDFLAGS += -no-undefined
endif
