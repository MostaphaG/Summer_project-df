from setuptools import setup, find_packages
 
classifiers = [
  'Development Status :: 5 - Production/Stable',
  'Intended Audience :: Education',
  'Intended Audience :: End Users/Desktop',
  'Operating System :: MacOS :: MacOS X',
  'Operating System :: Microsoft :: Windows',
  'License :: OSI Approved :: MIT License',
  'Programming Language :: Python :: 3'
]
 
setup(
  name='dformpy',
  version='1.0',
  description='Differential forms and exterior algebra visualizer',
  long_description=open('README.txt').read() + '\n\n' + open('CHANGELOG.txt').read(),
  long_description_content_type='text/plain',
  url='https://github.com/MostaphaG/Summer_project-df/tree/main/Library',  
  author='Moustafa Gharamti, Maciej Jarema, Samuel Kirwin-Jones',
  author_email='moustafa.gharamti@nottingham.ac.uk, macusjarema@gmail.com, 1999samkj@gmail.com',
  license='MIT', 
  classifiers=classifiers,
  keywords='Differential forms, Exterior algebra , Vector fields , Exterior derivative, Interior derivative, Hodge, Wedge , Curl , Gradiant , Divergence , Derivative',  
  packages=find_packages(),
  install_requires=[''] 
)
