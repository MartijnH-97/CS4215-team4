# Generator and scheduler of jobs

The sparkgen tool combines a job generator and scheduler supportiung multiple priotities (strict) with or without preeemption.
The tool uses a config file in json format.

## Compilation

Set the environment variables for go:

`setenv`

Compile via:

`go build sparkgen`

## Cmdline options

Supported options are:

- `-c` to specify the config file to use (default: config.json)
- `-d` to enable more debug output
- `-r` to record the stdout and stderr of the run jobs (will cerate a log directory with a pair of files for each scheduled job)

Example:

`./sparkgen -r -d -c conf.json.sample`

## Config

A sample configuration is conf.json.sample.

Th options are
- `jobClassParameters` a vector of structs to define jobs classes, specifying the system and hyper parameters of that class as well as its probabilites and priorites.
- `master` the addres of Spark's master node
- `lambda` the overall job arrival rate
- `runtime` the runtime of one experiment
- `preemptJobs` if to preemept jobs once a higher priority job is generaated or not