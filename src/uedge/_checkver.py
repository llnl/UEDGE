def _check_newer_uedge_ver():
    pkg = "uedge"

    try:
        import json
        import urllib.request

        try:
            import importlib.metadata

            thisver = importlib.metadata.version(pkg)
        except:
            import pkg_resources

            thisver = pkg_resources.get_distribution(pkg).version

        contents = urllib.request.urlopen("https://pypi.org/pypi/" + pkg + "/json").read()
        data = json.loads(contents.decode())
        thatver = data["info"]["version"]

        if thisver < thatver:
            return f"\nAn update to UEDGE v{thatver} is available via PyPi (pip)"
    except Exception as err:
        return "Error checking pypi version: {}".format(err)
