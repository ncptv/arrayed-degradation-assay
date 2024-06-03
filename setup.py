from setuptools import setup, find_packages

with open('requirements.txt') as f:
    required = f.read().splitlines()

setup(
    name='arrayed_degradation_assay',
    version='0.1.0',
    description='A package for arrayed degradation assay analysis.',
    url='https://github.com/ncptv/arrayed-degradation-assay',
    author='Inceptive Nucleics, Inc.',
    author_email='opensource@inceptive.team',
    license='Apache 2.0',
    packages=find_packages(),
    install_requires=required,
    python_requires='>=3.8',
)
