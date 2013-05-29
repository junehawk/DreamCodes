win3 = 0
win2 = 0
out3 = c(1:31)
out2 = c(1:31)
out3[26] = 0
out2[26] = 0
sum = 0
total = 0
min1 = 10000
max1 = -1
min2 = 10000
max2 = -1
mi1 = 0
mim1 = 0
mi2 = 0
mim2 = 0
for(i in 1:31) {
	if(i==26) {
		total = total + 1
		next
	}
	rmseFilename1 = paste("/Users/junehawk/Dream7/out3/Drug",i,"_ModelRMSE.txt", sep="")
	rmseloaded1 = read.table(rmseFilename1, check.names=F)
	rmses1 = t(rmseloaded1)
	orderedRMSE1 = rmses1[order(rmses1,decreasing=F)]
	rmseFilename2 = paste("/Users/junehawk/Dream7/out2/Drug",i,"_ModelRMSE.txt", sep="")
	rmseloaded2 = read.table(rmseFilename2, check.names=F)
	rmses2 = t(rmseloaded2)
	orderedRMSE2 = rmses2[order(rmses2,decreasing=F)]
	if(orderedRMSE1[1]<=orderedRMSE2[1]) {
		win3 = win3+1
	} else {
		win2 = win2 +1
	}
	out3[i] = orderedRMSE1[1]
	out2[i] = orderedRMSE2[1]

	for(j in 1:length(orderedRMSE1)) {
		if(orderedRMSE1[j]<10) {
			sum = sum+orderedRMSE1[j]
			total = total +1
		}
		if(orderedRMSE1[j]<min1) {
			min1 = orderedRMSE1[j]
		}
		if(orderedRMSE1[j]>max1 & orderedRMSE1[j]<10) {
			max1 = orderedRMSE1[j]
			mi1 = i
			mim1 = j
		}
	}
	for(j in 1:length(orderedRMSE2)) {
		if(orderedRMSE2[j]<min2) {
			min2 = orderedRMSE2[j]
		}
		if(orderedRMSE2[j]>max2 & orderedRMSE2[j]<10) {
			max2 = orderedRMSE2[j]
			mi2 = i
			mim2 = j
		}
	}
}
cat("Total model number is ", total , "\n")
cat("Average RMSE is ", sum/total)