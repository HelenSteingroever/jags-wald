ACLOCAL_AMFLAGS = -I m4
SUBDIRS = libltdl src win

.PHONY: win32-install win64-install installer

win32-install:
	make prefix=`pwd`/win/inst32
	make install prefix=`pwd`/win/inst32

win64-install:
	make prefix=`pwd`/win/inst64 
	make install prefix=`pwd`/win/inst64

installer:
	(cd win && $(MAKE) $@)

examples-zip:
	(zip -r examples examples)

FORCE:
