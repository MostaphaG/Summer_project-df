from setuptools import setup, find_packages
 
classifiers = [
  'Development Status :: 5 - Production/Stable',
  'Intended Audience :: Education',
  'Operating System :: Microsoft :: Windows :: Windows 10',
  'License :: OSI Approved :: MIT License',
  'Programming Language :: Python :: 3'
]
 
setup(
  name='formpy',
  version='1.0',
  description='Differential forms and Exterior Algebra visualizer',
  long_description=open('README.txt').read() + '\n\n' + open('CHANGELOG.txt').read(),
  url='https://github.com/MostaphaG/Summer_project-df/tree/main/FormPy_Library',  
  author='Moustafa Gharamti','Maciej Jarema', 'Samuel Kirwin-Jones'
  author_email='moustafa.gharamti@nottingham.ac.uk','macusjarema@gmail.com', '1999samkj@gmail.com'
  license='MIT', 
  classifiers=classifiers,
  keywords='calculator', 
  packages=find_packages(),
  install_requires=[''] 
)
