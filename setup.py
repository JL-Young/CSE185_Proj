from setuptools import setup, find_packages

setup(
    name='t_sKNEE',
    version='1.0.0',  # You can specify the version
    description='CSE1185 Final Project',
    author='Li, Venkatasubramani, Young',
    author_email='jal098@ucsd.edu, jvenkatasubramani@ucsd.edu, jlyoung@ucsd.edu',
    url='https://github.com/JL-Young/CSE185_Proj',
    packages=find_packages(),  # Automatically find packages in your project directory
    install_requires=[
        'numpy',
        'matplotlib',
        'scanpy',
        'leidenalg',
    ],
)
