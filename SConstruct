version = '0.4.0'

import os

AddOption('--prefix', dest='prefix', type='string', nargs=1, action='store',
          metavar='DIR', default='/usr',
          help='installation prefix')
AddOption('--libdir', dest='libdir', type='string', nargs=1, action='store',
          metavar='DIR', default='/usr/lib',
          help='library installation directory')
AddOption('--debug-build', action='store_true', default=False,
          help='build debug version')

cflags = '-D_GNU_SOURCE -Werror -Wall -Wextra -Wmissing-prototypes ' +\
         '-Winit-self -Wcast-align -Wpointer-arith ' +\
         '-Wno-unused-parameter -Wuninitialized -Wno-sign-compare'
libs = ['zxcvbn', 'm']
if GetOption('debug_build'):
    cflags += ' -g -O0 -fstack-protector-all ' +\
              '-fsanitize=undefined -fno-omit-frame-pointer -fsanitize=address'
    libs += ['ubsan', 'asan']
else:
    cflags += ' -O2'
cflags += ' ' + os.environ.get('CFLAGS', '')

env = Environment(PREFIX=GetOption('prefix'),
                  LIBDIR=GetOption('libdir'),
                  SHLIBVERSION=version)
if 'CC' in os.environ:
    env['CC'] = os.environ['CC']
if 'LIBPATH' in os.environ:
    env['LIBPATH'] = os.environ['LIBPATH'].split(':') + env.get('LIBPATH', [])

libzxcvbn = env.SharedLibrary('zxcvbn', 'zxcvbn.c',
                              CFLAGS=cflags)
zxcvbn_cli = env.Program('zxcvbn_cli', 'zxcvbn_cli.c',
                         LIBS=libs, LIBPATH=['.'] + env.get('LIBPATH', []),
                         CFLAGS=cflags)
Default([libzxcvbn, zxcvbn_cli])

env.Alias('install', [env.InstallVersionedLib('$LIBDIR', libzxcvbn),
                      env.Install('$PREFIX/include', 'zxcvbn.h')])
