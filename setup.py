from setuptools import setup
setup(
  name = 'myocardiumModel',         # How you named your package folder (MyLib)
  packages = ['myocardiumModel'],   # Chose the same as "name"
  version = '2.0.3',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'simplified myocardium motion model',   # Give a short description about your library
  author = 'Yu Zheng',                   # Type in your name
  author_email = 'jorry.zhengyu@gmail.com',      # Type in your E-Mail
  url = 'https://github.com/JorryZ/myocardiumModel',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/JorryZ/myocardiumModel/archive/v2.0.3.tar.gz',    # I explain this later on
  keywords = ['myocardium', 'motion', 'shell model'],   # Keywords that define your package best
  install_requires=['numpy','scipy','trimesh','meshplex'],
  classifiers=[
    'Development Status :: 5 - Production/Stable',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package    
    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',    
    'License :: OSI Approved :: MIT License',   # Again, pick a license    
    'Programming Language :: Python :: 3.6',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.7',
  ],
)
