import os
import operator
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict

small_np_to_time=defaultdict(list)
medium_np_to_time=defaultdict(list)
large_np_to_time=defaultdict(list)

def GetAxis(dictionary):
	xAxisContent=[]
	yAxisContent=[]
	for key in sorted(dictionary.iterkeys()):
		xAxisContent.append(key)
		yAxisContent.append(dictionary[key])
	return xAxisContent,yAxisContent



def GenerateLineChart(dict1,dict2,dict3,imageName1,imageName2,imageName3,xAxisLabel,yAxisLabel,fileName):
	x1,y1 = GetAxis(dict1)
	x2,y2 = GetAxis(dict2)
	x3,y3 = GetAxis(dict3)

	fig, ax = plt.subplots()
	plt.plot(x1,y1,marker="x")
	plt.plot(x2,y2,marker="x")
	plt.plot(x3,y3,marker="x")
	ax.set(xlabel=xAxisLabel, ylabel=yAxisLabel)	
	ax.grid()
	plt.legend([imageName1, imageName2,imageName3], loc=2)
	path="./data/graphs/"
	plt.savefig((path+fileName), format='eps', dpi=1000)
	#plt.show()


def initDictionary(dictionary, path):
	#print path
	# read all files and create a hashtable 
	# with key equals  to the number of threads
	# and values equals to list with the times
	for fileName in os.listdir(path):
		if(".txt" not in fileName):
			continue
		f = open(path+fileName,"r")
		lines =  f.readlines();
		for line in lines:
			if("total time" in line):
				time = float((line.split(" ")[3]).strip())
				np = int(line.split(" ")[7].strip())
				dictionary[np].append(time)

	# sort  time measurments
	for key in dictionary.keys():
		dictionary[key] = sorted(dictionary[key])
	#print dictionary
	# update hash table value for each key with median time instead of list of times
	for key in dictionary.keys():
		medianTime = (dictionary[key])[len(dictionary[key])/2]
		dictionary[key] =  medianTime

def InitSpeedUpDictionary(target,src):

	for key in src.keys():
		target[key] = src[1]/src[key]



def main():
	
	initDictionary(small_np_to_time,"./data/qsub/small/")
	initDictionary(medium_np_to_time,"./data/qsub/medium/")
	initDictionary(large_np_to_time,"./data/qsub/large/")
	#print small_np_to_time
	#print medium_np_to_time
	#print large_np_to_time
	GenerateLineChart(small_np_to_time,medium_np_to_time,large_np_to_time,"192x128","512x384","768x768","number of processes","total time in (sec)","time.eps")

	small_speedup = dict()
	medium_speedup = dict()
	large_speedup = dict()

	InitSpeedUpDictionary(small_speedup,small_np_to_time)
	InitSpeedUpDictionary(medium_speedup,medium_np_to_time)
	InitSpeedUpDictionary(large_speedup,large_np_to_time)

	print small_speedup
	print medium_speedup
	print large_speedup
	GenerateLineChart(small_speedup,medium_speedup,large_speedup,"192x128","512x384","768x768","number of processes","Speedup","speedup.eps")

if __name__ == "__main__":
	main()





