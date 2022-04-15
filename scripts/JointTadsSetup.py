from setuptools import setup
from distutils.extension import Extension
import numpy

setup(
    name="JointTads",
    ext_modules=[
        Extension("JointTadsLib",
                  ["JointTadsLib.pyx"],
                  include_dirs=[numpy.get_include()],
                  #define_macros=[("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION")],
                  language="c++")]
)