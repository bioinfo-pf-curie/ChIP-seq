# Nextflow profiles


If `-profile` is not specified at all the pipeline will be run locally and expects all software to be installed and available in the `PATH` environment variable. The following `-profile` are available:

## Define where the tools are available


###  `conda`

Build a new conda environment before running the pipeline. Use the option `--condaCacheDir` to change the default conda cache directory.

###  `multiconda`

Build a new conda environment for each process before running the pipeline. Use the option `--condaCacheDir` to change the default conda cache directory.

###  `path`

Use a global path for all tools. Use the option `--globalPath` to define the path the use.

###  `multipath`

Use the paths defined in configuration for each tool.

First, create a folder tree that looks like this 

```
├── path
│   ├── alpine
│   │   └── bin
│   ├── fastqc
│   │   └── bin
│   ├── helloWorld
│   │   └── bin
│   ├── rmarkdown
│   │   └── bin
│   └── trickySoftware
│       └── bin
```

The second level in the folder tree is the label that has been given in each Nextflow process. Install each tool in its corresponding `bin` directory.

Then, use the option `--globalPath` to define the path to the `path` folder as shown above in the folder tree.

###  `docker`

Use the Docker images for each process.

###  `singularity`

Use the Singularity images for each process. Use the option `--singularityImagePath` to specify where the images are available.

## Define where to launch the computation

By default, the pipline runs locally. If you want to run it on a computing cluster, use the profile below.

###  `cluster`

Run the workflow on the cluster, instead of locally.

## Test the pipeline

### `test`

Set up the test dataset
