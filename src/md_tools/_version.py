def get_version():
    # Get version from git tag
    try:
        from subprocess import check_output
        from shlex import split

        tag = check_output(split("git describe --tags --abbrev=0")).strip().decode("utf-8")
        return tag
    except Exception:  # noqa
        pass

    return "0.0.0"


VERSION = get_version()
