How to compile:

    >cd BACH/
    >make

How to run:

    All necessary information are stored in Configuration-file. See for example 
    xml/Configuration.xml
    Algorithms can be added or removed here, values of constatns can be changed
    without the need to recompile.
    To test the alignment algorithm, create a misaligned geometry first:

    >cd scripts
    >python MakeMisalign.py ../geom/Telescope_geom.xml
    >cd ..

    Now run the program:
    
    >bin/bach xml/Configuration.xml

    You can have a look at the output-histograms with

    >root -l out/Histograms.root
