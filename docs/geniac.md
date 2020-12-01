# Pipeline deployment with `geniac`

The `geniac` module allows the automatic generation of the [Nextflow](https://www.nextflow.io) configuration files (including profiles) and the container recipes and their creation.

## Prerequisite

* git (>= 2.0) [required]
* cmake (>= 3.0) [required]
* Nextflow (>= 20.01) [required]
* Singularity (>= 3.2) [optional]
* Docker (>= 18.0) [optional]

First, the user can defined the following environment variables:

```shell
#!/usr/bin/env bash
export WORK_DIR=/data/tmp/$USER/sandbox/nf-CRISPR
export SRC_DIR=$WORK_DIR/src
export INSTALL_DIR=$WORK_DIR/install
export BUILD_DIR=$WORK_DIR/build
```

## Initialization

Using the default behavior, only the installation folder where the pipeline will be deployed is required through the `CMAKE_INSTALL_PREFIX` option.

```shell
git clone --recursive https://gitlab.curie.fr/data-analysis/nf-CRISPR $SRC_DIR
mkdir -p $INSTALL_DIR $BUILD_DIR
cd $BUILD_DIR
cmake $SRC_DIR/geniac -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR
```

The `cmake` command can be used with the following options ;
* `ap_annotation_path`: Path to the annotation folder used by the pipeline. (`""`)
* `ap_install_docker_images`: Build the `docker` containers (`"OFF"`)
* `ap_install_docker_recipes`: Generate the `docker` recipes (`"OFF"`)
* `ap_install_singularity_images`: Build the `singularity` containers (`"OFF"`)
* `ap_install_singularity_recipes`: Gene the `singularity` recipes (`"OFF"`)
* `ap_nf_executor`: Cluster scheduler used by Nextflow (`"pbs"`)
* `ap_singularity_image_path`: Path to the `singularity` containers (`""`)
* `ap_use_singularity_image_link`: Use path to the `singularity` images (`"OFF"`)

## Example 

In the following example, the `singularity` images are available in the `SING_DIR` folder, and will thus be used during the  deployment of the pipeline 
using the options `ap_singularity_image_path` et `ap_use_singularity_image_link`.

```shell
cmake $SRC_DIR/geniac -DCMAKE_INSTALL_PREFIX=$INSTALL_DIR -Dap_singularity_image_path="${SING_DIR}" -Dap_use_singularity_image_link="ON"  
```

## Installation

Note that if the installation of `singularity` or `docker` images has been activated at the previous stage, the `make` command needs to be
run with the `root` rights.

```shell
cd $BUILD_DIR
make
make install
```

# Test

## Singularity


```shell
cd ${INSTALL_DIR}/pipeline  

nextflow run main.nf --singleEnd 'true' --genome 'hg38' --library 'GW-KO-Sabatini-Human-10' --samplePlan 'test/sample_plan.csv' -profile singularity
```

## Multiconda

The `multiconda` profile generates a `conda` environment for each tool listed in the `conf/geniac.conf` configuration.

```shell
cd ${INSTALL_DIR}/pipeline  

nextflow run main.nf --singleEnd 'true' --genome 'hg38' --library 'GW-KO-Sabatini-Human-10' --samplePlan 'test/sample_plan.csv' -profile multiconda
```
