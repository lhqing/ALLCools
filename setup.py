from setuptools import setup, find_packages

allcools_version = '0.1'

setup(
    name='allcools',
    version=allcools_version,
    author='Hanqing Liu',
    author_email='hanliu@salk.edu',
    packages=find_packages(),
    description='Tool kit for ALLC format and methylation analysis.',
    long_description=open('README.md').read(),
    include_package_data=True,
    install_requires=['pandas>=0.24', 'numpy', 'psutil'],
    entry_points={
        'console_scripts': ['allcools=ALLCools.__main__:main'],
    }
)

if __name__ == '__main__':
    f = open("ALLCools/__init__.py", 'w')
    f.write(f"__version__ = '{allcools_version}'\n")
    f.close()
