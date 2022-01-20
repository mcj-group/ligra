import os

Import('env', 'runtime')

env = env.Clone()
env.Append(CPPPATH = os.path.join(Dir('.').srcnode().abspath, 'ligra'))

SConscript(dirs=['apps'], exports = ['env', 'runtime'])
