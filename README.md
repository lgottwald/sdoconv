# sdoconv

Converts vensim mdl-files containing system dynamics models to gams
modelling language. Additionally the mdl files can be extended with vpd-files
to specify an objective function and with voc-files to add controls to
the model. A set of corresponding vpd-, voc- and mdl-files can be grouped
in a vop-files.

Specifications of the formats of the file formats and examples will be
added soon.

## Build/Install

To build the converter use the standard cmake work flow, i.e.
run the following commands from the source directory
```
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
```
And if you wish to install additionally run
```
make install
```

The converter only depends on the parser module libsdo which will be
automatically fetched during build.
The parser module depends on bison, flex and boost.