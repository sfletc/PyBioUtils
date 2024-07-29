from setuptools import setup, find_packages

setup(
    name='py_bio_utils',
    version='0.1.3',
    description='DNA sequence manipulation tools',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Stephen Fletcher',
    author_email='steveofbrisbane@gmail.com',
    license='MIT',
    url='https://github.com/sfletc/PyBioUtils',
    packages=find_packages(),  # Automatically discover and include all packages in the package directory
    install_requires=[
        'numpy>=1.20.0',
        'textwrap3>=0.9.2',
        'biopython'
    ],
    extras_require={
        'gzip_support': ['gzip-reader>=0.1.0']
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
    python_requires='>=3.6',
)
