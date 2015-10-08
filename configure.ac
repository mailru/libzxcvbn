AC_PREREQ([2.63])
AC_INIT([zxcvbn], [0.1], [k.kogdenko@corp.mail.ru])
AM_INIT_AUTOMAKE([1.9 foreign no-define])
AC_CONFIG_HEADERS([config.h])
: ${CFLAGS="-Wall -Wextra -Wmissing-prototypes -Wno-sign-compare -Wno-char-subscripts -Wno-unused-parameter -D_GNU_SOURCE"}
AC_PROG_CC
AM_PROG_LIBTOOL

AC_CHECK_HEADERS([sys/queue.h],, AC_MSG_ERROR("Couldn't find 'sys/queue.h'"))

AC_OUTPUT([Makefile])
