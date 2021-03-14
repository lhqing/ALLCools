from setuptools import setup, find_packages

setup(
    name='allcools',
    use_scm_version=True,
    setup_requires=['setuptools_scm'],
    author='Hanqing Liu',
    author_email='hanliu@salk.edu',
    description='Tool kit for single-cell methylome data analysis.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/lhqing/ALLCools',
    license='MIT',
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=find_packages(exclude=('docs', 'test')),
    include_package_data=True,
    install_requires=[],
    entry_points={
        'console_scripts': ['allcools=ALLCools.__main__:main'],
    }
)
