from pybuilder.core import Author, init, use_plugin

use_plugin("python.core")
use_plugin('exec')
use_plugin('copy_resources')

#use_plugin("python.unittest")
use_plugin("python.install_dependencies")
use_plugin("python.flake8")
use_plugin("python.coverage")
use_plugin("python.distutils")

use_plugin('pypi:pybuilder_header_plugin')
use_plugin('pypi:pybuilder_release_plugin')

name = "biocuration"
version = "0.5"

authors = [Author("kp14", "harvardsommer@gmx.net")]
license = 'Public Domain'
url = "https://bitbucket.org/kp14/biocuration"
description = __doc__
summary = "Tools for biocuration."

default_task = ["publish"]


@init
def set_properties(project):
    project.build_depends_on('coverage')
    project.set_property('coverage_break_build', False)

    project.build_depends_on('wheel')

    #project.set_property("dir_source_main_python", r"src\main\python\biocuration")
    project.set_property("dir_source_unittest_python", r"src/unittest/python")
    #project.set_property("unittest_module_glob", "test_*.py")
    #project.set_property("unittest_test_method_prefix", "test_")
    project.set_property("run_unit_tests_command",
    "py.test -v %s" % project.expand_path("$dir_source_unittest_python"))
    project.set_property("run_unit_tests_propagate_stdout", True)
    project.set_property("run_unit_tests_propagate_stderr", True)

    project.set_property('copy_resources_target', '$dir_dist')
    project.get_property('copy_resources_glob').append('LICENSE.txt')
    project.get_property('copy_resources_glob').append('setup.cfg')

    project.set_property('flake8_verbose_output', True)
    project.set_property('flake8_break_build', False)
    project.set_property('flake8_include_test_sources', True)

    project.set_property('pybuilder_header_plugin_break_build', False)
    project.set_property('pybuilder_header_plugin_expected_header', open('header.py').read())

    project.set_property('distutils_commands', 'bdist_wheel')
    project.set_property('distutils_classifiers', [
        "Development Status :: 3 - Alpha",
        'License :: Public Domain'])
