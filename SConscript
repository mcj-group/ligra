import os

Import('env', 'runtime')

env = env.Clone()
env.Append(CPPPATH = os.path.join(Dir('.').srcnode().abspath, 'ligra'))

# Just as with pbbs, this suite is riddled with warnings.
env.Append(CPPFLAGS = ['-Wno-all'])
# Several ligra apps (e.g. CF) can't handle gcc's -ftree-loop-vectorize pass,
# which comes with -O3. This manifests in segfaults. Fortunately clang doesn't
# accept that flag, so we don't need to disable it.
if not GetOption('clang'): env.Append(CPPFLAGS = ['-fno-tree-loop-vectorize'])

# [mcj] The applications won't build with the following defined.
# So I guess we're using 32-bit integers
#env.Append(CPPDEFINES = ['LONG', 'EDGELONG'])
# TODO decide about -DPD -DBYTE -DNIBBLE  -DBYTERLE

if runtime != "competition":
    env = env.Clone()
    env.Append(CPPDEFINES = ["CAN_USE_SWARM_API"])
    # Imitate the default Ligra parallelization strategy: use coarse grain inner
    # loops for vertices with degree less than one thousand.
    cgEnv = env.Clone()
    cgEnv.Append(CPPDEFINES = ['PLS_LIGRA_COARSE_GRAIN'])
    cgEnv['PROGSUFFIX'] = '_cg'
    cgEnv['OBJSUFFIX'] = '.cgo'
    fgEnv = env.Clone()
    fgEnv['PROGSUFFIX'] = '_fg'
    fgEnv['OBJSUFFIX'] = '.fgo'

    SConscript(dirs=['apps', 'apps/bucketing'],
               exports = ['fgEnv', 'cgEnv'])

else:
    serialEnv = env.Clone()
    parallelEnv = env.Clone()
    serialEnv['PROGPREFIX'] = 's'
    serialEnv['OBJSUFFIX'] = '.os'
    parallelEnv['PROGPREFIX'] = 'p'
    parallelEnv['OBJSUFFIX'] = '.op'

    hasCilk = False
    if not GetOption('no_exec'):
        conf = parallelEnv.Configure()
        hasCilk = conf.CheckLibWithHeader('cilkrts', 'cilk/cilk.h', 'c')
        parallelEnv = conf.Finish()
    if GetOption('no_exec') or hasCilk:
        parallelEnv.Append(CPPFLAGS = ['-fcilkplus'])
        parallelEnv.Append(CPPDEFINES = ['CILK'])

    SConscript(dirs=['apps', 'apps/bucketing'],
               name='SConscript-competition',
               exports = ['serialEnv', 'parallelEnv'])
