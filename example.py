def f():
    pass


code = f.__code__

print(code.co_name)
print(code.co_filename)
print(code.co_lnotab)
