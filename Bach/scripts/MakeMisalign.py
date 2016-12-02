from __future__ import print_function

import os
import sys

import xml.etree.ElementTree as ET
import random
from random import gauss

# Seed the generator
random.seed(42)


if len(sys.argv) != 2:
    print("usage: python MakeMisalign.py <geometry.xml>")
    sys.exit()

tree = ET.parse(sys.argv[1])
root = tree.getroot()
module_info = {}
thispath = os.path.abspath(__file__).replace('MakeMisalign.py', '')
outfile = open(thispath+"/../geom/Misalign_geom.xml", "write")
outfile.write("""<?xml version="1.0" encoding="utf-8"?>

<Geometry>
""")
for child in root:
    outfile.write("\t<Module")
    outfile.write(" name='"+child.attrib['name']+"'")
    outfile.write(" X='"+child.attrib['X']+"'")
    outfile.write(" Y='"+child.attrib['Y']+"'")
    outfile.write(" Z='"+child.attrib['Z']+"'")
    outfile.write(" RX='"+child.attrib['RX']+"'")
    outfile.write(" RY='"+child.attrib['RY']+"'")
    outfile.write(" RZ='"+child.attrib['RZ']+"'")
    if child.attrib['name'] == 'D09-W0108':
        outfile.write(" dX='"+str(float(child.attrib['dX'])+gauss(0, 0.0))+"'")
        outfile.write(" dY='"+str(float(child.attrib['dY'])+gauss(0, 0.0))+"'")
        outfile.write(" dZ='"+child.attrib['dZ']+"'")
        outfile.write(" dRX='"+child.attrib['dRX']+"'")
        outfile.write(" dRY='"+child.attrib['dRY']+"'")
        outfile.write(" dRZ='"+str(float(child.attrib['dRZ'])+gauss(0, 0.00))+"'")
    else:
        outfile.write(" dX='"+str(float(child.attrib['dX'])+gauss(0, 0.05))+"'")
        outfile.write(" dY='"+str(float(child.attrib['dY'])+gauss(0, 0.05))+"'")
        outfile.write(" dZ='"+child.attrib['dZ']+"'")
        outfile.write(" dRX='"+str(float(child.attrib['dRX'])+gauss(0, 0.00))+"'")
        outfile.write(" dRY='"+str(float(child.attrib['dRY'])+gauss(0, 0.00))+"'")
        outfile.write(" dRZ='"+str(float(child.attrib['dRZ'])+gauss(0, 0.01))+"'")
    outfile.write("/>\n")
outfile.write("</Geometry>")
outfile.close()
