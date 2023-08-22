from setuptools import setup

setup(
    name='py_bio_utils',
    version='0.1.0',
    description='DNA sequence manipulation tools',
    author='Stephen Fletcher',
    author_email='steveofbrisbane@gmail.com',
    license='MIT',
    url='https://github.com/sfletc/PyBioUtils',
    packages=['pi_bio_utils'],
    install_requires=[
        'numpy>=1.20.0',
        'textwrap3>=0.9.2', # note that the 'textwrap' module is part of Python's standard library, so we are using 'textwrap3' here
    ],
    extras_require={
        'gzip_support': ['gzip-reader>=0.1.0'] # optional, in case there is a separate package needed for gzip support
    },
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)