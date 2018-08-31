import xml.etree.ElementTree as ET
import os

###################################
# Goal is to create selection.h
###################################

###################################
# should become a parameter
###################################
ifilename='selectionDEC2016PU.xml'
template='template.h'
ofilename='moriondPU.h'

tree = ET.parse(ifilename)
#tree = ET.parse('selection2016.xml')
#tree = ET.parse('test.xml')
root = tree.getroot()


import SelectionMaker as sel
s = sel.Selection()

###################################
# Read Baseline
###################################
for child in root:
	if child.tag == "Baseline":
		s.SetBaseline(child.attrib['selection'])

###################################
# Add Variable
###################################
for var in root.iter("Variable"):
   #print var.attrib
   s.AddVariable(var.attrib)

###################################
# Add Bins
###################################
for bins in root.iter("Bin"):
   #print bins.attrib
   s.AddBin(bins.attrib)

###################################
# Add Region
###################################
for regions in root.iter("Region"):
   s.AddRegion(regions.attrib)

###################################
# Create selection region
###################################
template='template.h'
ofilename='moriondPU.h'
ofilename2="moriondPU.cc"
ofilename3="yieldStrMPU.h"
ofilename4="SRStrMPU.h"
ofilename5="TFStrMPU.h"
ofilename6="SystMPU.h"
command = "rm " + " " + ofilename +" " + ofilename2 +" " + ofilename3 +" " + ofilename4 +" " + ofilename5+" " + ofilename6
os.system(command)
command = "cp "+template+" "+ofilename
os.system(command)

s.CreateSelFunctions(ofilename)
s.AddSelection(ofilename2)
s.DumpAllRegionsVectors(ofilename3)
s.DumpSignalRegionsVectors(ofilename4)
s.DumpTFRegionsVectors(ofilename5)
s.DumpSignalSystRegionsVectors(ofilename6)


