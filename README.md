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

cd ~/MyAIDABach/build && cmake -DDD4hep_DIR=$DD4hep_DIR .. && make -j1 && cd .. && cd example/ && LD_LIBRARY_PATH=/afs/cern.ch/user/c/cburr/MyAIDABach/build/Bach/lib:$LD_LIBRARY_PATH geoPluginRun -volmgr -destroy -plugin Bach_main -config GenerateConfiguration.xml -input file:$HOME/MyAIDABach/example/Telescope.xml -deltas file:$HOME/MyAIDABach/example/repository.xml && LD_LIBRARY_PATH=/afs/cern.ch/user/c/cburr/MyAIDABach/build/Bach/lib:$LD_LIBRARY_PATH geoPluginRun -volmgr -destroy -plugin Bach_main -config Configuration.xml -input file:$HOME/MyAIDABach/example/Telescope.xml -deltas file:$HOME/MyAIDABach/example/repository-misalign.xml && cd .. && sha256sum example/geom/Geom_output_expected.xml example/output/Telescope_geom.xml && cat example/example_output.xml
