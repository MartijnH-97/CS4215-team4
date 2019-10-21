import os
import datetime

#lambdas = [0.00000001, 0.0000003]
#epochs = [3, 5, 7]
#nodes = [1, 2, 3, 4, 5]
#nReplications = 3

lambdas = [0.00000001]
epochs = [3]
nodes = [5]
nReplications = 1


now = datetime.datetime.now()
experiment = now.isoformat()

for replication in range(0, nReplications):
	for lambdaVal in lambdas:
		for epoch in epochs:
			for node in nodes:
				os.system("./src/sparkgen/sparkgen -d -c sample-conf.json -e {} -n {} -l {} -f {}".format(epoch, node, lambdaVal, experiment))

