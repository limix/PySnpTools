import os
import platform
import shutil
import sys
from distutils.command.clean import clean as Clean

import numpy
from setuptools import Extension, setup

# Version number
version = '0.3.12'


def readme():
    with open('README.md') as f:
        return f.read()


try:
    from Cython.Distutils import build_ext
except ImportError:
    use_cython = False
else:
    use_cython = True

#use_cython=False


class CleanCommand(Clean):
    description = "Remove build directories, and compiled files (including .pyc)"

    def run(self):
        Clean.run(self)
        if os.path.exists('build'):
            shutil.rmtree('build')
        for dirpath, dirnames, filenames in os.walk('.'):
            for filename in filenames:
                if (filename.endswith('.so') or filename.endswith('.pyd')
                        #or filename.find("wrap_plink_parser.cpp") != -1 # remove automatically generated source file
                        #or filename.find("wrap_matrix_subset.cpp") != -1 # remove automatically generated source file
                        or filename.endswith('.pyc')):
                    tmp_fn = os.path.join(dirpath, filename)
                    print("removing", tmp_fn)
                    os.unlink(tmp_fn)


# set up macro
if platform.system() == "Darwin":
    macros = [("__APPLE__", "1")]
elif "win" in platform.system().lower():
    macros = [("_WIN32", "1")]
else:
    macros = [("_UNIX", "1")]

#see http://stackoverflow.com/questions/4505747/how-should-i-structure-a-python-package-that-contains-cython-code
if use_cython:
    ext_modules = [
        Extension(
            name="pysnptools.snpreader.wrap_plink_parser",
            language="c++",
            sources=[
                "pysnptools/snpreader/wrap_plink_parser.pyx",
                "pysnptools/snpreader/CPlinkBedFile.cpp"
            ],
            include_dirs=[numpy.get_include()],
            define_macros=macros),
        Extension(
            name="pysnptools.snpreader.wrap_matrix_subset",
            language="c++",
            sources=[
                "pysnptools/snpreader/wrap_matrix_subset.pyx",
                "pysnptools/snpreader/MatrixSubset.cpp"
            ],
            include_dirs=[numpy.get_include()],
            define_macros=macros)
    ]
    cmdclass = {'build_ext': build_ext, 'clean': CleanCommand}
else:
    ext_modules = [
        Extension(
            name="pysnptools.snpreader.wrap_plink_parser",
            language="c++",
            sources=[
                "pysnptools/snpreader/wrap_plink_parser.cpp",
                "pysnptools/snpreader/CPlinkBedFile.cpp"
            ],
            include_dirs=[numpy.get_include()],
            define_macros=macros),
        Extension(
            name="pysnptools.snpreader.wrap_matrix_subset",
            language="c++",
            sources=[
                "pysnptools/snpreader/wrap_matrix_subset.cpp",
                "pysnptools/snpreader/MatrixSubset.cpp"
            ],
            include_dirs=[numpy.get_include()],
            define_macros=macros)
    ]
    cmdclass = {}


class CleanCommand(Clean):
    description = "Remove build directories, and compiled files (including .pyc)"

    def run(self):
        Clean.run(self)
        if os.path.exists('build'):
            shutil.rmtree('build')
        for dirpath, dirnames, filenames in os.walk('.'):
            for filename in filenames:
                if (filename.endswith('.so') or filename.endswith('.pyd') or
                        filename.find("wrap_plink_parser.cpp") !=
                        -1  # remove automatically generated source file
                        or filename.find("wrap_matrix_subset.cpp") !=
                        -1  # remove automatically generated source file
                        or filename.endswith('.pyc')):
                    tmp_fn = os.path.join(dirpath, filename)
                    print("removing", tmp_fn)
                    os.unlink(tmp_fn)


def setup_package():
    src_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    old_path = os.getcwd()
    os.chdir(src_path)
    sys.path.insert(0, src_path)

    needs_pytest = {'pytest', 'test', 'ptr'}.intersection(sys.argv)
    pytest_runner = ['pytest-runner>=2.9'] if needs_pytest else []
    setup_requires = ["cython", "numpy"] + pytest_runner
    tests_require = ['pytest']

    setup(
        name='pysnptools',
        version=version,
        description='PySnpTools',
        long_description=readme(),
        keywords='gwas bioinformatics sets intervals ranges regions',
        url=
        "http://research.microsoft.com/en-us/um/redmond/projects/mscompbio/",
        author='MSR',
        author_email='fastlmm@microsoft.com',
        license='Apache 2.0',
        packages=[
            "pysnptools/snpreader", "pysnptools/kernelreader",
            "pysnptools/pstreader", "pysnptools/standardizer",
            "pysnptools/kernelstandardizer", "pysnptools/util", "pysnptools"
        ],
        package_data={
            "pysnptools": [
                "test/datasets/all_chr.maf0.001.N300.bed",
                "test/datasets/all_chr.maf0.001.N300.bim",
                "test/datasets/all_chr.maf0.001.N300.fam",
                "test/datasets/phenSynthFrom22.23.N300.randcidorder.txt",
                "tests/datasets/all_chr.maf0.001.covariates.N300.txt"
            ]
        },
        install_requires=['scipy>=0.15.1', 'numpy>=1.9.2', 'pandas>=0.16.2'],
        setup_requires=setup_requires,
        tests_require=tests_require,
        cmdclass=cmdclass,
        ext_modules=ext_modules)

    try:
        setup(**metadata)
    finally:
        del sys.path[0]
        os.chdir(old_path)


if __name__ == '__main__':
    setup_package()
