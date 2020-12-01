# Installation

This documentation has been adapted from the nf-core guidelines
(see https://nf-co.re/usage/installation for details).

To start using this pipeline, follow the steps below:

1. [Install Nextflow](#install-nextflow)
2. [Install the pipeline](#install-the-pipeline)
3. [Geniac](#geniac)
4. [Pipeline configuration](#pipeline-configuration)
    * [Resources](#resources)
    * [Software dependencies](#software-dependencies)
5. [Cluster usage](#cluster-usage)
6. [Reference genomes](#reference-genomes)

## Install Nextflow

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
nextflow run main.nf
```

If you would like to make changes to the pipeline, it's better to fork the repository on your github account and then clone the pipeline from your personal repository. 
Once cloned, you can run the pipeline directly as above.

## Geniac

This current version of the pipeline is compatible with the `geniac` utility for automatic production deployment.
See the [`docs/geniac.md`](geniac.md) page for details and [geniac](https://geniac.readthedocs.io).

## Pipeline configuration

The default [`conf/base.config`](../conf/base.config) and  [`conf/process.config`](../conf/process.config) configuration are loaded by the pipeline.

They define some default parameters in particular the minimal computing resource requirements to run the processes. They are suitable for running the pipeline on a simple local server.
Note that a few variables related to [software dependencies](#software-dependencies) can be changed in this configuration file.

Be aware of two important points about this default configuration:

1. The default profile uses the `local` executor
    * All jobs run in the login session. If you use a simple server, this may be fine. 
	If you use a computing cluster, take care of not running all the jobs on the submission node.
    * See the [nextflow docs](https://www.nextflow.io/docs/latest/executor.html) for information about running with other hardware backends.
	Most job scheduler systems are natively supported.
2. Nextflow will expect all software to be installed and available on the `PATH` environment variable.
    * It's expected to use an additional config profile for [conda](https://docs.conda.io), [singularity](https://sylabs.io/guides/3.6/user-guide/) and [Docker](https://www.docker.com/) (See below).

### Resources

The required resources for each process are defined in the configuration [`conf/process.config`](../conf/process.config).
These resources are usually defined using the `label` directive associated with the number of RAM/CPUs available for each job.
This file can be edited if one wants to change the resources allocated to a specific process.

### Software dependencies

Software dependencies can be managed in several ways. We provide here an overview and we recommend the user to read carefully the [Profiles](profiles.md) section to properly use the different possibilities.

#### Paths

By default, Nextflow expects all the tools to be installed and available in the `PATH` environment variable.
This path can be set using the `-profile path` combined with the `--globalPath /my/path` option specifying where the tools are installed.
In addition, the `-profile multipath` is available in order to specify the `PATH` for each tool, instead of a global one.

#### Conda

If you're not able to use [singularity](https://sylabs.io/guides/3.6/user-guide/) _or_ [Docker](https://www.docker.com/), you can instead use [conda](https://docs.conda.io) to manage the software requirements.
This is slower and less reproducible than the above, but is still better than having to install all requirements yourself!
The pipeline ships with a conda environment file and Nextflow has a built-in support for this.

To use it first ensure that you have [conda](https://docs.conda.io) installed (we recommend to use [miniconda](https://conda.io/miniconda.html)), then follow the same pattern as above and use the flag `-profile conda`
Note that in this case, the conda environment will be created in the `$HOME/conda-cache-nextflow` folder by default. This folder can be changed using the `--condaCacheDir` option.

In addition to a general conda environment, this pipeline also comes with a `-profile multiconda` setting. In this case, a conda environment per tool (note that each process is assigned a tool with the label directive) will be created.
This configuration is useful if different processes require different tool versions (leading to potential conda conflicts thus making impossible to use the `-profile conda` option).

#### Singularity

Using [singularity](https://sylabs.io/guides/3.6/user-guide/) is in general a great idea to manage environment and ensure reproducibility.
The process is very similar: running the pipeline with the option `-profile singularity` tells Nextflow to enable singularity to run the pipeline. 
Containers with all of the software requirements can be automatically generated using the `recipes` information.
Once available, the user can specified where to look for the images using the option `--singularityImagePath /my/path/to/singularity/containers`.

#### Docker

A generic configuration profile is available with `-profile docker`.
In this case, the pipeline will look for `Docker` containers as defined in the [`conf/docker.config`](conf/docker.config).

### Cluster usage

By default, we set up a `cluster` profile to execute the pipeline on a computational cluster.
Please, edit the `cluster.config` file to set up your own cluster configuration.

### Reference genomes

See [`docs/referenceGenomes.md`](referenceGenomes.md).
