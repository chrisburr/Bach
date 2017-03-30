# Bach

Software for aligning telescope like detectors

## Building

### Setup

```bash
export DD4hep_DIR=$HOME/DD4hep-2017
export BACH_DIR=$HOME/MyAIDABach

# Source DD4hep
cd ~/DD4hep-2017
source bin/thisdd4hep.sh
```

### Compiling
```bash
mkdir -p $BACH_DIR/build
cd $BACH_DIR/build
cmake -D DD4hep_DIR=$DD4hep_DIR ..
make
```

## Running the example
```bash
cd $BACH_DIR/example/
export LD_LIBRARY_PATH=$BACH_DIR/build/Bach/lib:$LD_LIBRARY_PATH

geoPluginRun -volmgr -destroy -plugin Bach_main -config GenerateToy.xml -input file:$BACH_DIR/example/xml/Telescope.xml -deltas file:$BACH_DIR/example/xml/repository-nominal.xml

geoPluginRun -volmgr -destroy -plugin Bach_main -config AlignToyData.xml -input file:$BACH_DIR/example/xml/Telescope.xml -deltas file:$BACH_DIR/example/xml/repository-misaligned.xml
```

![AIDA-2020 Acknowledgement](http://aida2020.web.cern.ch/sites/aida2020.web.cern.ch/files/AIDA-2020%20acknowledgement%20presentations_0.png)
