from setuptools import setup
setup(
  name = 'TranscriptClean', # the name of your package
  packages = ['TranscriptClean'], # same as above
  version = '2.1', # version number
  license='MIT', # license type
  description = 'TranscriptClean corrects common long-read RNA-seq artifacts', # short description
  author = 'Dana Wyman', # your name
  author_email = 'fairlie.reese@gmail.com', # your email
  url = 'https://github.com/mortazavilab/TranscriptClean/', # url to your git repo
  download_url = 'https://github.com/mortazavilab/TranscriptClean/archive/refs/tags/v2.1.tar.gz', # link to the tar.gz file associated with this release
  keywords = ['transcriptclean', 'transcription', 'isoform'], #
  install_requires=[ # these can also include >, <, == to enforce version compatibility
      'pybedtools>0.8',
      'pyfaidx>0.5',
      'pytest>6',
      'pyranges==0.0.127'
      ],
  classifiers=[ # choose from here: https://pypi.org/classifiers/
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Science/Research ',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3.7',
  ],
  entry_points={
    "console_scripts": [
        'transcriptclean=TranscriptClean.TranscriptClean:main',
        'transcriptclean_get_sjs=TranscriptClean.accessory_scripts.get_SJs_from_gtf:main']}
)
