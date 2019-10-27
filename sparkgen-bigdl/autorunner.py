import os
import datetime

lambdas = [0.002, 0.006]
epochs = [1, 5]
nodes = [1, 3]
nReplications = 1

#lambdas = [0.00277777777]
#epochs = [2]
#nodes = [3, 5]
#nReplications = 1

now = datetime.datetime.now()
experiment = now.isoformat()

for replication in range(0, nReplications):
	for lambdaVal in lambdas:
		for epoch in epochs:
			for node in nodes:
				os.system("./src/sparkgen/sparkgen -d -c conf.json -e {} -n {} -l {} -f {}".format(epoch, node, lambdaVal, experiment))
		os.system("y | ./remove-cluster.sh")
		os.system("./create-cluster.sh")


