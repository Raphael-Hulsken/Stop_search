Here I give some explanation for the program that I used/create :

-raphael.cc : It the main file that I used, it is used to do all the selection of events, to apply the differents region of signal. The main of the code is in the ../sonicScrewdriver/interface/BabyScrewdriver.h file. If you want to have information about the cut and the region go to the ./Selection/moriond.h file. And if you want to have information about the variables go to the ../common/Reader_CommonFormat_CommonBabies.h file. The output are the plot in the plotsTest/ folder and the cards with the information about the number of events per region. The cards could be used with other code such as resultsPlot.cc combined_comparaison.cc or Zbi_result.cc .

-resulPlot.cc : this code use the card that are the output of the raphael.cc program. It is used to plot the number of events per region or the Zbi per region. By rewritting the code you should be able to plot the both in the same plots. You should take care that this code need, at the moment I write this Readme file, 3 card of data that are in the folder /table_of_data and 2 file where you have the name of the regions, in the folder name_of_region, that are obtained with the raphael.cc code.  You should take care of what you want for table of data and put the right name of region with it otherwise it will not work.

-combined_comparaison.cc : same utilisation as resulPlot.cc but here it write as output the number of events per regions. I used it in order to compare the events in R0 and R1 and to be sure that we have in total the same number of events. For this code you need one card and one file with the name of the region.

-Zbi_result.cc : the use is the same as the 2 files above. It give at output the Zbi for every region. 

-compare_bkg.cc : This code is create in order to compare the data in the code with the data that is in the 2016 paper that I wrote in some file. It give as an output the ratio between the paper and the simulation

-createLimiPlot.cc : This code will plot the limit plots and the ratio cross section mesured and cross section excluded. It need root plots that are writen by the combine software. You need to check cards1/combineCommands.sh to see how it works and the explaination are given in the twiki (ask others for that twiki). For know it needs 4 root file in input but it could be easily changed. You could use the file that are in the dat_limit/ or dat_limit2/ folder in order to understand how it work but I advise you to write your own root file quickly =)

-normalisation.cc : at the beggining I used this code to plot the same plot as the one schedule with raphael.cc but only for the signal in the stack way in order to see something to compare the different scenario. But now I use this code to see the efficiancy of the cut. It plots the ratio S/sqrt(B) with cut divided by S/sqrt(B) with cut.

-find_data.cpp : This code is just some test about the combine sofware. It was used to compare the value for some couple of masses.


-the TestSignificance/ folder : in this folder I make some code in order to understand the differences between the combine results and the Zbi results for differents configuration (such as different incertainties and different value for the signal) in function of the background. The C code are written in order to make the execution and to compute and the python code are here to plot the results with the value that you obtained with the C code. For each C code there should be a python code that plot the results

-the cards.../ folder : there are create by the combine maker of raphael.cc and it is used to compute the combine but it is only for one set of region per folder

- the dat_limit/ and dat_limit2/ : they are the folder where I put the root file for the limit. In the dat_limit it is the one with the old incertainties and in the second there are the new one

- the table_of_data/ folder : it is the one were I put the data file in but take care when you compile the raphael.cc code the file don't go in this folder but in the main folder

- the name_of_region/ folder : it is the one were I put the file with the name of the regions in but take care when you compile the raphael.cc code the file don't go in this folder but in the main folder

I you have questions about some code I write feel free to ask. =)
