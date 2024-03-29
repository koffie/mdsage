## -*- encoding: utf-8 -*-
import os
import re
from setuptools import setup
from codecs import open  # To open the README file with proper encoding
from distutils.command import build as build_module

if "SAGE" in os.environ:
    sage = os.environ["SAGE"]
else:
    sage = "sage"


# Obtain the different Sage versions
def get_all_version_names(
    mirror_url, idx=None, distribution="Ubuntu_14.04-x86_64", broken=[]
):
    import urllib2

    if idx is None:
        idx = 0
    else:
        idx = int(idx)
    site = urllib2.urlopen(mirror_url).read()
    ans = re.findall("(sage-([0-9]*(?:\.[0-9]*)*)-%s.tar.bz2)" % distribution, site)
    all_version_names = []
    for fname, ver in ans:
        if fname not in all_version_names:
            if not any(fname.startswith(i) for i in broken):
                all_version_names.append(fname)
    return all_version_names[idx]


# Get information from separate files (README, VERSION)
def readfile(filename):
    with open(filename, encoding="utf-8") as f:
        return f.read()


# Check the right Sage version
class build(build_module.build):
    def run(self):
        try:
            from sagemath.check_version import check_version

            check_version(sage_required_version)
        except ImportError:
            print("WARNING sagemath not installed, not checking version")
        build_module.build.run(self)


if __name__ == "__main__":
    cmdclass = {"build": build}

    # set the following to true for using cython to build the pyx files
    build_cython = False

    if not build_cython:
        ext_modules = []

    else:
        cython_build_ext = True
        # The next block is needed if there are cython files
        from setuptools import Extension

        try:
            from Cython.Build import cythonize
            import Cython.Compiler.Options
        except ImportError:
            cython_build_ext = False
        from sage.env import sage_include_directories

        ext_modules = [
            Extension(
                "mdsage.one_cython_file",
                sources=[os.path.join("mdsage", "one_cython_file.pyx")],
                include_dirs=sage_include_directories(),
            )
        ]
        if cython_build_ext:
            # Cython modules
            from Cython.Distutils import build_ext

            cmdclass["build_ext"] = build_ext
        else:
            # C files generated by making sdist
            ext_modules = [
                Extension(
                    "mdsage.one_cython_file",
                    sources=[os.path.join("mdsage", "one_cython_file.c")],
                    include_dirs=sage_include_directories(),
                )
            ]
            ext_modules = cythonize(ext_modules)

    # Specify the required Sage version
    sage_required_version = ">=7.4"
    REQUIREMENTS = [
        i.strip() for i in open("requirements.txt").readlines() if i[0] != "#"
    ]
    CLASSIFIERS = [
        i.strip() for i in open("classifiers.txt").readlines() if i[0] != "#"
    ]

    setup(
        name="mdsage",
        version=readfile(
            "VERSION"
        ),  # the VERSION file is shared with the documentation
        description="A repository for all the sage functions I've written. See: www.sagemath.org.",
        long_description=readfile(
            "README.rst"
        ),  # get the long description from the README
        url="https://github.com/koffie/mdsage",
        author="Maarten Derickx",
        author_email="maarten@mderickx.nl",  # choose a main contact email
        license="GNU General Public License v3",
        classifiers=CLASSIFIERS,
        keywords="SageMath, Number Theory, Modular Curves",
        install_requires=REQUIREMENTS,  # This ensures that Sage is installed
        packages=["mdsage"],
        ext_modules=ext_modules,  # This line is only needed if there are cython or c files present
        package_data={"": ["data_files/*"]},
        include_package_data=True,
        cmdclass=cmdclass,  # adding a special setup command for tests
    )
