import nox


@nox.session(python="3.9")
def tests(session):
    session.install(".[test]")
    hr_deps = ["pygeoogc", "pygeoutils"]
    for p in hr_deps:
        session.install(f"git+https://github.com/cheginit/{p}.git")
    session.run("pytest")
    session.run("coverage", "report")
    session.run("coverage", "html")
