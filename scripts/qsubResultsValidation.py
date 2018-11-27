import os
import filecmp
#B145772
#compares all images with the expected image
processNumber = [1,2,3,4,8,16,32]
expectedImage="./data/output/imagenew192x128_expectedOutput.pgm"
for i in processNumber:
	targetImg ="./data/qsub/small/imagenew192x128_"+str(i)+".pgm"
	res =  filecmp.cmp(targetImg,expectedImage);
	if(res == False):
		print "File compare failed for image: "+targetImg
		exit(-1);
	
print "File compares passed successfully!"

# valates that the average pixel is the same for every log file.
avgNum =-1
currAvgNum=-1
for i in processNumber:
	logFile ="./data/qsub/small/"+str(i)+".txt"
	f = open(logFile,"r")
	lines =  f.readlines();
	for line in lines:
		if("avg number" in line):
			currAvgNum =   line.split("=")[1].split(",")[0]
			if avgNum == -1:
				avgNum = currAvgNum
			else:
				assert(avgNum == currAvgNum)

print "Average pixel validation passed successfully!"