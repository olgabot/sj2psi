import sys

from setuptools import setup
from setuptools import find_packages


if sys.version_info[0] == 3:
    LONG_DESCRIPTION = open('README', encoding='utf-8').read()
else:
    LONG_DESCRIPTION = open('README').read()

setup(
    name='sj2psi',
    version='0.0.1',
    author='Olga B. Botvinnik',
    author_email='olga.botvinnik@gmail.com',
    packages=find_packages(),
    license='MIT License',
    url='http://olgabot.github.io/prettyplotlib',
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.3',
        'Topic :: Scientific/Engineering'
    ],
    description='Convert RNA-STAR SJ.out.tab files to 5-prime and 3-prime '
                '"percent spliced in" ("psi") scores.',
    long_description=LONG_DESCRIPTION,
    install_requires=['pandas>=0.13.1']
)
