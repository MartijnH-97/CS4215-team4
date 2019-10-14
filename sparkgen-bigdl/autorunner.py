import os

lambdas = [0.00000001, 0.0000003]
epochs = [3, 5, 7]
nodes = [1, 2, 3, 4, 5]
nReplications = 3

for replication in range(0, nReplications):
	for lambdaVal in lambdas:
		for epoch in epochs:
			for node in nodes:
				os.system("./src/sparkgen/sparkgen -d -c sample-conf.json -e {} -n {} -l {}".format(epoch, node, lambdaVal))

