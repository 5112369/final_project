from setuptools import setup, find_packages

setup(
    name='strobealign',
    version='0.1.0',
    author='liu zhengtao',
    author_email='5112369.@sjtu.edu.cn',
    description='A strobe alignment algorithm for bioinformatics',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/5112369/final_project',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'biopython',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)