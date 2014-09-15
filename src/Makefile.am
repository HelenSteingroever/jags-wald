SUBDIRS = distributions functions

jagsmod_LTLIBRARIES = bernoulli.la

bernoulli_la_SOURCES = bernoulli.cc

bernoulli_la_CPPFLAGS = -I$(includedir)/JAGS

bernoulli_la_LIBADD = distributions/bernoullidist.la
bernoulli_la_LIBADD += functions/bernoullifunc.la

if WINDOWS
bernoulli_la_LIBADD += -ljags-$(JAGS_MAJOR) -ljrmath-0
else
bernoulli_la_LIBADD += -ljags -ljrmath
endif

bernoulli_la_LDFLAGS = -module -avoid-version
if WINDOWS
bernoulli_la_LDFLAGS += -no-undefined
endif
