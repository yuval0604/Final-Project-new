from setuptools import setup, Extension

module = Extension('symnmf_capi', sources=['symnmfmodule.c', 'symnmf.c'])

setup(
    name='symnmf_capi',
    version='1.0',
    description='Symmetric NMF C API',
    ext_modules=[module],
)
