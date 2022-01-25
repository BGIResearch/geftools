geftools : Tools for manipulating GEFs
==============

For a full documentation, see [geftools GitHub page](https://bgiresearch.github.io/geftools/).


## Installation
To install stereopy from source, you need:
- HDF5 1.12.1 or newer with development headers
- OpenCV 4.5.4 or newer with development headers
- A C compiler

You can specify build options for geftools as environment variables when you build it from source:
```shell
git clone https://github.com/BGIResearch/geftools.git
cd geftools
HDF5_ROOT=/path/to/hdf5 OpenCV_DIR=/path/to/opencv cmake .
make
make install
```

The supported build options are:
- HDF5_ROOT: To specify where to find HDF5.
- OpenCV_DIR: To specify where to find OpenCV.
- CMAKE_INSTALL_PREFIX: Install directory used by make install.
- GEFTOOLS_BUILD_DOC: Option to build documentation.


## LIST OF COMMANDS
```text
Command: bgef          Generate common bin GEF(.bgef) according to gem file or bin1 GEF
         cgef          Generate cell bin GEF(.cgef) according to common bin GEF and mask file
         view          View GEF
```


## COMMANDS AND OPTIONS
### geftools bgef [OPTION...]
```text
  -i, --input-file FILE   input gene expression matrix file(.gem/.gem.gz) or bin1 bGEF file [request]
  -o, --output-file FILE  output bin GEF file (.bgef) [request]
  -b, --bin-size STR      Set bin size by the comma-separated list [request] (default: 1,10,20,50,100,200,500)
  -r, --region STR        Restrict to a rectangular region. The region is represented by the comma-separated list of
                          two vertex coordinates (minX,maxX,minY,maxY) (default: "")
  -t, --threads INT       number of threads (default: 8)
  -v, --verbose           Verbose output
      --help              Print help
```


### geftools cgef [OPTION...]
Generate cell bin GEF (.cgef) according to common bin GEF (.bgef) file and mask file
```text
  -i, --input-file FILE   input bin GEF file [request]
  -m, --mask-file FILE    input mask file [request]
  -o, --output-file FILE  output cell bin GEF file (.cgef) [request]
  -b, --block FILE        Pre block size (default: 256,256)
  -v, --verbose           Verbose output
      --help              Print help
```


### geftools view [OPTION...]
Show the contents of cell bin GEF
```text
-i, --input-file FILE     Input cell bin GEF file [request]
-o, --output-gem FILE     Output cell bin gem file (default: stdout)
-m, --output-mask FILE    Output border of polygons to mask format file (default: "")
-r, --region STR          Restrict to a rectangular region. The region is represented by the comma-separated list of
two vertex coordinates (minX,minY,maxX,maxY) (default: "")
-g, --genes [^]STR        Comma separated list of genes to include (or exclude with "^" prefix)
-G, --genes-file [^]FILE  File of genes to include (or exclude with "^" prefix))
--force-genes         Only warn about unknown subset genes
-v, --verbose             Verbose output
--help                Print help
```


## Coding Style Guide
See [here](docs/coding_style_guide.md).


