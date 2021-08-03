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

## Dashboard

The `moPLOT` package comes with a [`shiny`](https://shiny.rstudio.com)-based dashboard integrating the main visualizations and a wide variety of benchmark suites, allowing for an interactive exploration of their multi-objective landscapes. It can be started using the `runDashboard` command:

```r
library(moPLOT)
runDashboard()
```

## Citation

If you use our package to visualize your multi-objective problems using
* our gradient-field heatmaps (see our [EMO 2017 paper](http://link.springer.com/chapter/10.1007/978-3-319-54157-0_23) for details),
* our multi-objective PLOTs (see our [PPSN 2020 paper](https://link.springer.com/chapter/10.1007%2F978-3-030-58115-2_11) ([arXiv](https://arxiv.org/abs/2006.11547)) for details), or
* the `moPLOT` dashboard and the 3D-visualizations (see our [EMO 2021 paper](https://link.springer.com/chapter/10.1007/978-3-030-72062-9_50) ([arXiv](https://arxiv.org/abs/2011.14395)) for details),

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
  author="Sch{\"a}permeier, Lennart and Grimme, Christian and Kerschke, Pascal",
  editor="B{\"a}ck, Thomas and Preuss, Mike and Deutz, Andr{\'e} and Wang, Hao and Doerr, Carola and Emmerich, Michael and Trautmann, Heike",
  title="One PLOT to Show Them All: Visualization of Efficient Sets in Multi-objective Landscapes",
  booktitle="Parallel Problem Solving from Nature -- PPSN XVI",
  year="2020",
  publisher="Springer International Publishing",
  address="Cham",
  pages="154--167",
  isbn="978-3-030-58115-2"
}

@InProceedings{schaepermeier2021dashboard,
  author="Sch{\"a}permeier, Lennart and Grimme, Christian and Kerschke, Pascal",
  editor="Ishibuchi, Hisao and Zhang, Qingfu and Cheng, Ran and Li, Ke and Li, Hui and Wang, Handing and Zhou, Aimin",
  title="To Boldly Show What No One Has Seen Before: A Dashboard for Visualizing Multi-objective Landscapes",
  booktitle="Evolutionary Multi-Criterion Optimization",
  year="2021",
  publisher="Springer International Publishing",
  address="Cham",
  pages="632--644",
  isbn="978-3-030-72062-9"
}
```

## Contact

If you have any suggestions or ideas, or run into problems while running the code, please
use the [issue tracker](https://github.com/kerschke/moPLOT/issues) or send me an e-mail (<kerschke@uni-muenster.de>).
