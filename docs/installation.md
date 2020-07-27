# Installation

This documentation has been modified from the nf-core guidelines
(see https://nf-co.re/usage/installation for details).

To start using this pipeline, follow the steps below:

1. [Install Nextflow](#1-install-nextflow)
2. [Install the pipeline](#2-install-the-pipeline)
3. [Pipeline configuration](#3-pipeline-configuration)
    * [Define your own configuration](#31-define-our-own-configuration)
    * [Cluster usage](#32-cluster-usage)
    * [Software deps: Docker and Singularity](#33-software-deps-singularity)
    * [Software deps: Bioconda](#34-software-deps-conda)
4. [Reference genomes](#4-reference-genomes)

## NextFlow
Nextflow runs on most POSIX systems (Linux, Mac OSX etc). It can be installed by running the following commands:

```bash
# Make sure that Java v8+ is installed:
java -version

# Install Nextflow
curl -fsSL get.nextflow.io | bash

# Add Nextflow binary to your PATH:
mv nextflow ~/bin/
# OR system-wide installation:
# sudo mv nextflow /usr/local/bin
```

See [nextflow.io](https://www.nextflow.io/) for further instructions on how to install and configure Nextflow.

## Install the pipeline

You just need to download/clone the source code and transfer the pipeline files manually:

```bash
wget https://mypipeline/archive/master.zip
mkdir -p ~/mypipelines
unzip master.zip -d ~/mypipelines/
cd ~/mypipelines
nextflow run ~/mypipelines/mypipeline-master
```

If you would like to make changes to the pipeline, it's best to make a fork on your account and then clone the files. 
Once cloned you can run the pipeline directly as above.

## Geniac

This current version of the pipeline is compatible with the `geniac` utility for automatic production deployment.  
See the [`docs/geniac.md`](geniac.md) page for details.

## Pipeline configuration

By default, the pipeline loads a basic server configuration [`conf/base.config`](../conf/base.config) 
and a [`standard`](../conf/standard.config) profile.
This uses a number of sensible defaults for process requirements and is suitable for running
on a simple local server.

Note that a few variables related to software dependencies can be changed in this configuration file.

Be warned of two important points about this default configuration:

1. The default profile uses the `local` executor
    * All jobs are run in the login session. If you're using a simple server, this may be fine. 
	If you're using a compute cluster, take care of not running all jobs on the head node.
    * See the [nextflow docs](https://www.nextflow.io/docs/latest/executor.html) for information about running with other hardware backends.
	Most job scheduler systems are natively supported.
2. Nextflow will expect all software to be installed and available on the `PATH`
    * It's expected to use an additional config profile for docker, singularity or conda support. See below.

### Ressources

The required resources for each process are defined in the configuration [`conf/process.config`](../conf/process.config).
These resources are usually defined using `label` associated with the number of RAM/CPUs available for each job.
This file can be edited if one want to change the resources allocated to a specific process.

### Software dependencies

#### Paths

As mentioned below, by default, Nextflow expects all software to be installed and available on the `PATH`.  
This path can be set using the `-profile path` and the `--globalPath PATH` option specifying where the tools are installed.  
In addition, the `-profile multipath` is available in order to specify the `PATH` for each tool, instead of a global one.

#### Conda

If you're not able to use Docker _or_ Singularity, you can instead use conda to manage the software requirements.  
This is slower and less reproducible than the above, but is still better than having to install all requirements yourself!
The pipeline ships with a conda environment file and nextflow has built-in support for this.

To use it first ensure that you have conda installed (we recommend [miniconda](https://conda.io/miniconda.html)), then follow the same pattern as above and use the flag `-profile conda`
Note that in this case, the environment will be created in the `cache/work` folder. This folder can be changed using the ``--condaCacheDir` option.

In addition to a general conda environment, this pipeline also comes with a `-profile multiconda` setting. In this case, a conda environment per process will be generated.  
This configuration is useful if different processes require different tool versions (leading to potential conda conflicts).

#### Singularity

Using [Singularity](http://singularity.lbl.gov/) is in general a great idea to manage environment and ensure reproducibility.  
The process is very similar: running the pipeline with the option `-profile singularity` tells Nextflow to enable singularity for this run.  
Images containing all of the software requirements can be automatically generated using the `recipes` information.  
Once available, the user can specified where to look for the images using the option `--singularityImagePath PATH`

#### Docker

A generic configuration profile is available with `-profile docker`.   
In this case, the pipeline will look for `Docker` images as defined in the [`conf/docker.config`](conf/docker.config)

### Cluster usage

By default, we set up a `cluster` profile to execute the pipeline on a computational cluster.  
Please, edit the `cluster.config` file to set up your own cluster configuration.

### Reference genomes

See [`docs/reference_genomes.md`](reference_genomes.md)
