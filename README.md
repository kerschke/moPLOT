# moPLOT: Visualizing Multi-Objective Optimization Problems

## Description

Although many (real-world) optimization problems actually contain multiple objectives and
hence actually are multi-objective optimization problems, people often prefer to handle
them as single-objective problems (as those are easier to grasp).

In an attempt to make multi-objective problems more conceivable, *moPLOT* provides the
means to visualize bi- and tri-objective problems (i.e., problems with two or three
objectives). This is achieved by plotting the landscapes of their multi-objective gradients
using various visualization tools.


## Installation

Currently, *moPLOT* is only available within this development version, however, we of
course also plan to submit it to CRAN in the near future.

In the mean time, feel free to use the development version of this package:

```r
install.packages("devtools")
devtools::install_github("kerschke/moPLOT")
```


## Quickstart

Examples for the different visualizations using the current version of this package can be found in the [examples](/examples) folder.


## Citation

If you use our package to visualize your multi-objective problems using
* our multi-objective PLOTs (see our [PPSN 2020 paper](https://link.springer.com/chapter/10.1007%2F978-3-030-58115-2_11) ([arXiv](https://arxiv.org/abs/2006.11547)) for details), or
* our gradient-field heatmaps (see our [EMO 2017 paper](http://link.springer.com/chapter/10.1007/978-3-319-54157-0_23) for details),

please cite the respective publications.

```
@inproceedings{KerschkeGrimme2017Expedition,
  author    = {Pascal Kerschke and Christian Grimme},
  title     = {{An Expedition to Multimodal Multi-Objective Optimization Landscapes}},
  booktitle = {{Proceedings of the 9$^{th}$ International Conference on Evolutionary Multi-Criterion Optimization (EMO)}},
  pages     = {329~--~343},
  series    = {{Lecture Notes in Computer Science (LNCS)}},
  volume    = {11411},
  editor    = {Heike Trautmann and G{\"u}nter Rudolph and Kathrin Klamroth and Oliver Sch{\"u}tze and Margaret Wiecek and Yaochu Jin and Christian Grimme},
  year      = {2017},
  publisher = {Springer},
  address   = {M{\"u}nster, Germany},
  isbn      = {978-3-319-54157-0},
  doi       = {10.1007/978-3-319-54157-0_23},
  url       = {http://link.springer.com/chapter/10.1007/978-3-319-54157-0_23}
}

@InProceedings{schaepermeier2020plot,
  Title     = {{One PLOT to Show Them All: Visualization of Efficient Sets in Multi-Objective Landscapes}},
  Author    = {Schaepermeier, Lennart and Grimme, Christian and Kerschke, Pascal},
  Booktitle = {{Proceedings of the 16th International Conference on Parallel Problem Solving from Nature (PPSN XVI)}},
  Year      = {2020},
  Month     = {September},
  Pages     = {},
  Location  = {Leiden, The Netherlands},
  Publisher = {{Springer}},
  Note      = {accepted}
}
```


## Contact

If you have any suggestions or ideas, or run into problems while running the code, please
use the [issue tracker](https://github.com/kerschke/moPLOT/issues) or send me an e-mail (<kerschke@uni-muenster.de>).
