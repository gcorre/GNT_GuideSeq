[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.16759028.svg)](https://doi.org/10.5281/zenodo.16759028) [![Static Badge](https://img.shields.io/badge/License-custom-blue?style=flat)](./LICENSE.md)

------------------------------------------------------------------------

This repository consists of a snakemake workflow dedicated to the analysis of libraries derived from the GUIDE-Seq protocol, including iGUIDE-Seq, GUIDE-Seq2, Tag-Seq, OliTag-Seq or any other libraries based on the original protocol by Tsai et al. 2015.

The key features of the pipeline are :

-   High versatility with analysis of multiplexed libraries including samples from different organisms, gRNA and PAM sequences, PCR orientations ... in a single run.

-   Consideration of bulges and mismatches between the gRNA and gDNA alignments.

-   Consideration of multi-hit reads (assignment to random position among best hits).

-   Prediction of Off-target sites using SWOFFinder

-   Annotation to gene and oncogene (user defined lists)

-   Reporting with QC and OT quantification.

    ------------------------------------------------------------------------

# License

This work is distributed under an **Academic Non-Commercial License**. Commercial use is **not permitted** without explicit written consent from the authors. See the [LICENSE](./LICENSE.md) file for details.

------------------------------------------------------------------------

# Citation

If you use this pipeline in your research, please cite:

> Repository : <https://github.com/gcorre/GNT_GuideSeq>
>
> DOI: 10.5281/zenodo.16759028
>
> Publication :

# Documentation

For installation and usage instructions, please read the complete [documentation](https://raw.githack.com/gcorre/GNT_GuideSeq/refs/heads/master/docs/documentation.html).
