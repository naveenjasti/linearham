env = Environment(
    CPPPATH = ['eigen', 'src'])

# Doubles compilation time.
#env.Append(CCFLAGS='-O3 -msse2')

env.Program(
    target='test',
    source=Split(' src/core.cpp src/germline.cpp src/linalg.cpp src/smooshable.cpp src/test.cpp'))

