# Bach

Software for aligning telescope like detectors

### Building
```bash
mkdir build
cmake -DDD4hep_DIR=$DD4hep_DIR ..
make
```

### Running the example
```bash
cd example
LD_LIBRARY_PATH=/afs/cern.ch/user/c/cburr/MyAIDABach/build/Bach/lib:$LD_LIBRARY_PATH geoPluginRun -volmgr -destroy -plugin Bach_main Configuration.xml
```

![AIDA-2020 Acknowledgement](http://aida2020.web.cern.ch/sites/aida2020.web.cern.ch/files/AIDA-2020%20acknowledgement%20presentations_0.png)
