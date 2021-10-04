import setuptools

setuptools.setup(
    name="interlayer_energies_demo",
    version="0.1.0",
    url="https://github.com/qmc-hamm/interlayer_energies_demo",
    license='MIT License',
    author="Mick Krongchon, Lucas Wagner",
    author_email="lkwagner@illinois.edu",
    description="QMC data for the interaction of graphene bilayers",
    long_description=open('README.md').read(),
    packages=['interlayer_energies_demo'],
    install_requires=[],
    classifiers=[
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'License :: OSI Approved :: MIT License',
    ],
    include_package_data=True,
    package_data={'interlayer_energies_demo': ['data/*.csv']},
)
