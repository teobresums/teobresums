# Run a test
import foo_module

# memory test routines
# http://fa.bianp.net/blog/2013/different-ways-to-get-memory-consumption-or-lessons-learned-from-memory_profiler/
import os
import psutil
import resource

def memory_usage_psutil():
    process = psutil.Process(os.getpid())
    mem = process.get_memory_info()[0] / float(2 ** 20)
    return mem

def memory_usage_resource():
    rusage_denom = 1024.
    mem = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / rusage_denom
    return mem

if 1:
    for i in range(7):
        a = foo_module.foo(i, 42.)
        print(a)
else:
    # memory test
    for i in range(100):
        a = foo_module.foo(10000000, 42.)
        print(memory_usage_psutil(),memory_usage_resource())
        #del a

