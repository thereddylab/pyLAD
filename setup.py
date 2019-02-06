#!/usr/bin/env python
"""pyLAD: Python-based DamID Lamina-Associated-Domain analysis"""

import os
import sys
import site
import subprocess
import setuptools


DOCLINES = __doc__.split("\n")

if sys.version_info[:2] < (2, 6) or (3, 0) <= sys.version_info[0:2] < (3, 4):
    raise RuntimeError("Python version 2.6, 2.7 or >= 3.4 required.")

if sys.version_info[0] < 3:
    import __builtin__ as builtins
else:
    import builtins
    from importlib import reload

CLASSIFIERS = """\
Intended Audience :: Science/Research
Intended Audience :: Developers
Programming Language :: Fortran
Programming Language :: Python
Programming Language :: Python :: 3
Topic :: Software Development
Topic :: Scientific/Engineering
Operating System :: Microsoft :: Windows
Operating System :: POSIX
Operating System :: Unix
Operating System :: MacOS

"""

MAJOR = 1
MINOR = 0
MICRO = 0
ISRELEASED = True
VERSION = '%d.%d.%d' % (MAJOR, MINOR, MICRO)


# Return the git revision as a string
def git_version():
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'Fortran'
        env['LANG'] = 'Fortran'
        env['LC_ALL'] = 'Fortran'
        out, err = subprocess.Popen(cmd, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE, env=env
                                    ).communicate()
        if err.startswith(b'fatal'):
            raise OSError()
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


# BEFORE importing distutils, remove MANIFEST. distutils doesn't properly
# update it when the contents of directories change.
if os.path.exists('MANIFEST'):
    os.remove('MANIFEST')


def get_version_info():
    # Adding the git rev number needs to be done inside
    # write_version_py(), otherwise the import of scipy.version messes
    # up the build under Python 3.
    FULLVERSION = VERSION
    if os.path.exists('.git'):
        GIT_REVISION = git_version()
    elif os.path.exists('pyLAD/version.py'):
        # must be a source distribution, use existing version file
        # load it as a separate module to not load scipy/__init__.py
        import imp
        version = imp.load_source('pyLAD.version',
                                  'pyLAD/version.py')
        GIT_REVISION = version.git_revision
    else:
        GIT_REVISION = "Unknown"

    if not ISRELEASED:
        FULLVERSION += '.dev-' + GIT_REVISION[:7]

    return FULLVERSION, GIT_REVISION


def write_version_py(filename='pyLAD/version.py'):
    cnt = """
# THIS FILE IS GENERATED FROM pyLAD SETUP.PY
short_version = '%(version)s'
version = '%(version)s'
full_version = '%(full_version)s'
git_revision = '%(git_revision)s'
release = %(isrelease)s

if not release:
    version = full_version
"""
    FULLVERSION, GIT_REVISION = get_version_info()

    a = open(filename, 'w')
    try:
        a.write(cnt % {'version': VERSION,
                       'full_version': FULLVERSION,
                       'git_revision': GIT_REVISION,
                       'isrelease': str(ISRELEASED)})
    finally:
        a.close()


try:
    from sphinx.setup_command import BuildDoc
    HAVE_SPHINX = True
except ImportError:
    HAVE_SPHINX = False

if HAVE_SPHINX:
    class pyLADBuildDoc(BuildDoc):
        """Run in-place build before Sphinx doc build"""
        def run(self):
            ret = subprocess.call([sys.executable, sys.argv[0],
                                   'build_ext', '-i'])
            if ret != 0:
                raise RuntimeError("Building pyLAD failed!")
            BuildDoc.run(self)


def configuration(parent_package='', top_path=None):
    DISTNAME = 'pyLAD'
    from numpy.distutils.misc_util import Configuration
    config = Configuration(
        DISTNAME, parent_package, top_path,
        #description='Python-based LAD analysis tool',
        #maintainer='Michael Sauria',
        #maintainer_email='mike.sauria@gmail.com',
        #url='https://github.com/thereddylab/pyLAD',
        scripts=['bin/LADetector'],
        )
    config.add_data_files('pyLAD/__init__.py')
    config.add_extension("dnacopy", sources=[
                         "pyLAD/cbststats.f",
                         "pyLAD/changepoints.f",
                         "pyLAD/getbdry.f",
                         "pyLAD/lchoose.f",
                         "pyLAD/phyper.f",
                         "pyLAD/pnorm.f",
                         "pyLAD/segmentp.f",
                         "pyLAD/smoothCNA.f",
                         "pyLAD/sort.f",
                         "pyLAD/tailprobs.f"])
    config.add_extension("string_trimming", sources=[
                         "pyLAD/string_trimming.f"])
    config.get_version('pyLAD/version.py')

    return config


def setup_package():

    # Rewrite the version file every time
    write_version_py()
    FULLVERSION, GIT_REVISION = get_version_info()

    if HAVE_SPHINX:
        cmdclass = {'build_sphinx': ScipyBuildDoc}
    else:
        cmdclass = {}

    metadata = dict(
        maintainer="Michael Sauria",
        maintainer_email="mike.sauria@gmail.com",
        description=DOCLINES[0],
        # long_description = "\n".join(DOCLINES[2:]),
        url='https://github.com/thereddylab/pyLAD',
        download_url='https://github.com/thereddylab/pyLAD/tarball/%s' % (
            FULLVERSION),
        license='BSD',
        cmdclass=cmdclass,
        package_dir={'': './'},
        install_requires=['numpy', 'scipy', 'pysam'],
        scripts = ['bin/LADetector'],
        classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
        platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
        test_suite='nose.collector',
        zip_safe=False,
    )

    if (len(sys.argv) >= 2
        and ('--help' in sys.argv[1:]
        or sys.argv[1] in ('--help-commands', 'egg_info', '--version',
                           'clean'))):
        # For these actions, NumPy is not required.
        #
        # They are required to succeed without Numpy for example when
        # pip is used to install Scipy when Numpy is not yet present in
        # the system.
        try:
            from setuptools import setup
        except ImportError:
            from distutils.core import setup

        metadata['version'] = FULLVERSION
    else:
        #################
        # This section is very hacky, but numpy.distutils.core.setup
        # doesn't seem to respect the 'install_requires' argument
        # so we use a series of dummy installs to require each package
        # and get them installed.
        #
        # We need to use numpy.distutils.cor.setup in order to get
        # the fortran code to be automatically compiled.
        #################
        # Make sure numpy gets installed if it isn't already
        try:
            import numpy
        except ImportError:
            setuptools.setup(**dict(
                install_requires=['numpy'],
                name='pyLAD',
                version=FULLVERSION,
                classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
                platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
                ))
            # Reload 'site' to ensure that sys.path has numpy in it
            reload(site)
        # Make sure scipy gets installed if it isn't already
        try:
            import scipy
        except:
            setuptools.setup(**dict(
                install_requires=['scipy'],
                name='pyLAD',
                version=FULLVERSION,
                classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
                platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
                ))
            # Reload 'site' to ensure that sys.path has scipy in it
            reload(site)
        # Make sure pysam gets installed if it isn't already
        try:
            import pysam
        except:
            setuptools.setup(**dict(
                install_requires=['pysam'],
                name='pyLAD',
                version=FULLVERSION,
                classifiers=[_f for _f in CLASSIFIERS.split('\n') if _f],
                platforms=["Windows", "Linux", "Solaris", "Mac OS-X", "Unix"],
                ))
            # Reload 'site' to ensure that sys.path has pysam in it
            reload(site)

        from numpy.distutils.core import setup

        cwd = os.path.abspath(os.path.dirname(__file__))
        metadata['packages'] = setuptools.find_packages(
            exclude=['examples', 'test', 'ez_setup.py'], where='./')
        metadata['configuration'] = configuration

    setup(**metadata)


if __name__ == '__main__':
    setup_package()
