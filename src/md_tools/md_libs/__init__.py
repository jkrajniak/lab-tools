
try:
    from . import bonds  # noqa
except (ImportError, ValueError):
    print(" No bonds module")
