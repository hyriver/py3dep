import nox


@nox.session(python="3.10")
def tests(session):
    session.install(".[test,dem]")
    hr_deps = ["pygeoogc", "pygeoutils"]
    for p in hr_deps:
        session.install(f"git+https://github.com/cheginit/{p}.git", "--use-feature=in-tree-build")
    session.run("pytest")
    session.run("coverage", "report")
    session.run("coverage", "html")
