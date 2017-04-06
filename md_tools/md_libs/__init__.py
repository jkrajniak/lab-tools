from . import files_io
try:
    from . import bonds
except (ImportError, ValueError):
    print(' No bonds module')
