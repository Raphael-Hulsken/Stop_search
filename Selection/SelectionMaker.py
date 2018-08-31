
signs_ = ['>','>=','<',"<=",'==']

class Cut:
  
  def __init__(self):
     self.label = ""
     self.babytuple_name = ""
     self.cut_value = 0
     self.sign = ""


  def load(label, babytupleName, selection):
     self.label = label
     self.babytuple_name = babytupleName



  #def CreateSelFunction():


  #def AddSelection():



class Selection:


   def __init__(self):
	self.variables = []
	self.bins = []
	self.regions = []
	self.baseline = ""

   def AddVariable(self, variable_from_xml):
	self.variables.append(variable_from_xml)

   def AddBin(self, bin_from_xml):
	# change the selection with the following correspondance
	# and -> &&
	# !>  -> <
	bin = bin_from_xml
	bin['selection'] = bin['selection'].replace('and','&&', 100) # max nof occurence = 100
	bin['selection'] = bin['selection'].replace('or','||', 100) # max nof occurence = 100
	bin['selection'] = bin['selection'].replace('!>','<', 100) # max nof occurence = 100
   	self.bins.append(bin)
   
   def SetBaseline(self, baseline):
	self.baseline = baseline
	self.baseline = self.baseline.replace('and','&&', 100) # max nof occurence = 100
	self.baseline = self.baseline.replace('or','||', 100) # max nof occurence = 100
	self.baseline = self.baseline.replace('!>','<', 100) # max nof occurence = 100



   def AddRegion(self, region_from_xml):
	region = region_from_xml
	region['selection'] = region['selection'].replace('and','&&', 100) # max nof occurence = 100
	region['selection'] = region['selection'].replace('or','||', 100) # max nof occurence = 100
	region['selection'] = region['selection'].replace('!>','<', 100) # max nof occurence = 100
	self.regions.append(region)


   def CreateSelFunctions(self, ofilename):
     with open(ofilename, "a") as myfile:


	

	## replace MET by its babytuple name
	METbabytuple = "MET"
	for var in self.variables:
	  if var['name'] == "MET":
	    METbabytuple = "myEvent."+var['babyTuple']
	####################################
	
	
	####################################
	# Create Baseline function
	####################################
	baseSel = self.baseline
	for var in self.variables:
		baseSel = baseSel.replace(var['name'],"myEvent."+var['babyTuple'])
	dump = "bool baseline(){ return "+baseSel+";}\n"
	myfile.write(dump)
	print dump

	####################################
	# Create Region function
	####################################
	for region in self.regions:
	       selection = region['selection']
	       #replace the variable label by their babytuple name
	       for var in self.variables:
	         selection = selection.replace(var['name'],"myEvent."+var['babyTuple'])
	   
	       dump = "bool "+region['name']+"() { return ("+selection+" && baseline()  ); }\n" 
	       myfile.write(dump)
	       print dump


	################################
	## loop over the regions
	################################
	for region in self.regions:

	  ################################
	  ## loop over the bins
	  ################################
	  for bin in self.bins:
	     METbins = bin['METbins'].split(",")
	   
	     ################################
	     ## loop over the MET bins
	     ################################
	     for i  in range(len(METbins)-1):
	       #if i!=len(METbins)-1:
	       name=region['name']+bin['name']+METbins[i]
	       selection = region['name']+"() && "+bin['selection']+" && "+METbabytuple+">="+METbins[i]
	       if METbins[i+1]!="inf" :
	         selection+=" && "+METbabytuple+"<"+METbins[i+1]  
		 name+="lessMETless"+METbins[i+1]
	       else: 
	         name+="lessMETlessInf"
	       
	       #replace the variable label by their babytuple name
	       for var in self.variables:
	         selection = selection.replace(var['name'],"myEvent."+var['babyTuple'])
	      	 
	       dump="bool "+name+"() { return ("+selection+");}\n"
	       myfile.write(dump)
	       print dump

   def AddSelection(self,ofilename):
     with open(ofilename, "a") as myfile:
	################################
	## loop over the regions
	################################
	for region in self.regions:

	  ################################
	  ## loop over the bins
	  ################################
	  for bin in self.bins:
	     METbins = bin['METbins'].split(",")
	   
	     ################################
	     ## loop over the MET bins
	     ################################
	     for i  in range(len(METbins)-1):
	       #if i!=len(METbins)-1:
	       name=region['name']+bin['name']+METbins[i]
	       if METbins[i+1]!="inf" :
		 name+="lessMETless"+METbins[i+1]
	       else: 
	         name+="lessMETlessInf"
	      	 
	       dump="AddRegion(\""+name+"\",\""+name+"\",&"+name+");\n" #and all the regions for the up/down , do this more clever with a list
	       myfile.write(dump)
	       print dump
               #good ones systs = ["LSFdown", "LSFup", "BTlightDown", "BTlightUp","BTheavyDown", "BTheavyUp", "PUdown", "PUup", "PDFdown", "PDFup", "alphaSdown", "alphaSup", "Q2down", "Q2up", "ISRnjetsDown", "ISRnjetsUp" ]
               #systs = ["LSFdown", "LSFup", "BTlightDown", "BTlightUp","BTheavyDown", "BTheavyUp", "PUdown", "PUup", "PDFdown", "PDFup", "alphaSdown", "alphaSup", "Q2down", "Q2up", "topPtModelingdown", "topPtModelingup" ]
               #systs = ["PDFdown", "PDFup", "alphaSdown", "alphaSup", "Q2down", "Q2up" ]
               systs = []
               for syst in systs:
	           dump="AddRegion(\""+name+syst+"\",\""+name+syst+"\",&"+name+");\n"
	           myfile.write(dump)


   def DumpAllRegionsVectors(self,ofilename):
     with open(ofilename, "a") as myfile:
	################################
	## loop over the regions
	################################
	dump="vector<string> yield = { "
	myfile.write(dump)
	for region in self.regions:

	  ################################
	  ## loop over the bins
	  ################################
	  for bin in self.bins:
	     METbins = bin['METbins'].split(",")
	   
	     ################################
	     ## loop over the MET bins
	     ################################
	     for i  in range(len(METbins)-1):
	       #if i!=len(METbins)-1:
	       name=region['name']+bin['name']+METbins[i]
	       if METbins[i+1]!="inf" :
		 name+="lessMETless"+METbins[i+1]
	       else: 
	         name+="lessMETlessInf"
	      	 
	       dump= "\""+ name + "\""
	       myfile.write(dump)
               #if bin not len(self.bins):
	       dump= " , "
	       myfile.write(dump)

        dump = " };\n"
        myfile.write(dump)

   def DumpSignalRegionsVectors(self,ofilename):
     with open(ofilename, "a") as myfile:
	################################
	## loop over the regions
	################################
	dump="vector<string> signalReg = { "
	myfile.write(dump)
	for region in self.regions:

	  ################################
	  ## loop over the bins
	  ################################
	  for bin in self.bins:
	     METbins = bin['METbins'].split(",")
	   
	     ################################
	     ## loop over the MET bins
	     ################################
	     for i  in range(len(METbins)-1):
	       #if i!=len(METbins)-1:
	       name=region['name']+bin['name']+METbins[i]
	       if METbins[i+1]!="inf" :
		 name+="lessMETless"+METbins[i+1]
	       else: 
	         name+="lessMETlessInf"
	      
               substring = "SR1l"	 
               if substring in name:
	         dump= "\"" +name +"\""
	         myfile.write(dump)
	         dump= " , "
	         myfile.write(dump)

        dump = " };\n"
        myfile.write(dump)


   def DumpTFRegionsVectors(self,ofilename):
     with open(ofilename, "a") as myfile:
	################################
	## loop over the regions
	################################
	dump="vector<string> tfreg = { "
	myfile.write(dump)
	for region in self.regions:

	  ################################
	  ## loop over the bins
	  ################################
	  for bin in self.bins:
	     METbins = bin['METbins'].split(",")
	   
	     ################################
	     ## loop over the MET bins
	     ################################
	     for i  in range(len(METbins)-1):
	       #if i!=len(METbins)-1:
	       name=region['name']+bin['name']+METbins[i]
	       if METbins[i+1]!="inf" :
		 name+="lessMETless"+METbins[i+1]
	       else: 
	         name+="lessMETlessInf"
	      
               substring = "SR1l"	 
               if substring in name:
                 name2 = name.replace("SR1l", "")
	         dump= "\"" + name2 + "\"" 
	         myfile.write(dump)
	         dump= " , "
	         myfile.write(dump)

        dump = " };\n"
        myfile.write(dump)

   def DumpSignalSystRegionsVectors(self,ofilename):
     with open(ofilename, "a") as myfile:
	################################
	## loop over the regions
	################################
	dump="vector<string> yield = { "
	myfile.write(dump)
	for region in self.regions:

	  ################################
	  ## loop over the bins
	  ################################
	  for bin in self.bins:
	     METbins = bin['METbins'].split(",")
	   
	     ################################
	     ## loop over the MET bins
	     ################################
	     for i  in range(len(METbins)-1):
	       #if i!=len(METbins)-1:
	       name=region['name']+bin['name']+METbins[i]
	       if METbins[i+1]!="inf" :
		 name+="lessMETless"+METbins[i+1]
	       else: 
	         name+="lessMETlessInf"
	      
               substring = "SR1l"	 
               if substring in name:
	         dump= "\"" +name +"\""
	         myfile.write(dump)
	         dump= " , "
	         myfile.write(dump)
                 # good ones systs = ["LSFdown", "LSFup", "BTlightDown", "BTlightUp","BTheavyDown", "BTheavyUp", "PUdown", "PUup", "PDFdown", "PDFup", "alphaSdown", "alphaSup", "Q2down", "Q2up", "ISRnjetsDown", "ISRnjetsUp"  ]
                 #systs = ["LSFdown", "LSFup", "BTlightDown", "BTlightUp","BTheavyDown", "BTheavyUp", "PUdown", "PUup", "PDFdown", "PDFup", "alphaSdown", "alphaSup", "Q2down", "Q2up", "topPtModelingdown", "topPtModelingup" ]
                 #systs = ["PDFdown", "PDFup", "alphaSdown", "alphaSup", "Q2down", "Q2up" ]
                 systs = []
                 for syst in systs:
	             dump= "\"" +name+syst+"\""
	             myfile.write(dump)
	             dump= " , "
	             myfile.write(dump)



        dump = " };\n"
        myfile.write(dump)
