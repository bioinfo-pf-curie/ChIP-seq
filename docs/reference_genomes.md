# Reference Genomes Configuration

The pipeline needs a reference genome for alignment and annotation.
All annotation data and paths must be defined/modified in the `conf/genomes.conf` files.

These paths can be supplied on the command line at run time (see the [usage docs](../usage.md)),
but for convenience it's often better to save these paths in a nextflow config file.
See below for instructions on how to do this.

## Adding paths to a config file
Specifying long paths every time you run the pipeline is a pain.
To make this easier, the pipeline comes configured to understand reference genome keywords which correspond
to preconfigured paths, meaning that you can just specify `--genome ID` when running the pipeline.

Note that this genome key must be specified in the config file `genomes.conf`.

To use this system, add paths to your config file using the following template:

```nextflow
params {
  genomes {
    'YOUR-ID' {
      fasta  = '<PATH TO FASTA FILE>/genome.fa'
    }
    'OTHER-GENOME' {
      // [..]
    }
  }
  // Optional - default genome. Ignored if --genome 'OTHER-GENOME' specified on command line
  genome = 'YOUR-ID'
}
```

You can add as many genomes as you like as long as they have unique IDs.



