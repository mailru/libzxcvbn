version = '0.4.0'

import os

common_cflags = '-g -O2 -D_GNU_SOURCE -Wall -Wextra -Wmissing-prototypes ' + \
                '-Wno-sign-compare -Wno-char-subscripts -Wno-unused-parameter'

AddOption('--prefix', dest='prefix', type='string', nargs=1, action='store',
          metavar='DIR', default='/usr',
          help='installation prefix')
AddOption('--libdir', dest='libdir', type='string', nargs=1, action='store',
          metavar='DIR', default='/usr/lib',
          help='library installation directory')

env = Environment(PREFIX=GetOption('prefix'),
                  LIBDIR=GetOption('libdir'),
                  SHLIBVERSION=version,
                  CFLAGS=os.environ.get('CFLAGS',''))

libzxcvbn = env.SharedLibrary('zxcvbn', 'zxcvbn.c',
                              CFLAGS=common_cflags + ' $CFLAGS')
zxcvbn_cli = env.Program('zxcvbn_cli', 'zxcvbn_cli.c',
                         LIBS=['zxcvbn', 'm'], LIBPATH='.',
                         CFLAGS=common_cflags + ' $CFLAGS')
Default([libzxcvbn, zxcvbn_cli])

env.Alias('install', [env.InstallVersionedLib('$LIBDIR', libzxcvbn),
                      env.Install('$PREFIX/include', 'zxcvbn.h')])
