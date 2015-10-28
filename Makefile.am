noinst_bin_PROGRAMS=zxcvbn_cli
zxcvbn_cli_SOURCES=zxcvbn_cli.c
zxcvbn_cli_LDADD=libzxcvbn.la

lib_LTLIBRARIES=libzxcvbn.la
libzxcvbn_la_SOURCES=zxcvbn.c
libzxcvbn_la_LDFLAGS=-lm -version-number $(subst .,:,$(PACKAGE_VERSION))

include_HEADERS=zxcvbn.h

noinst_bindir=$(pwd)
