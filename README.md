## Add-GP-Bandits
This is a Matlab implementation of Additive Gaussian Process Optimisation and 
Bandits. For more details please read our paper (below).

### Installation & Getting Started
- Just add all relevant paths to your Matlab workspace
- The demos directory has two examples on how to use this library.
- demo.m: Easy set up and uses all default configurations.
- demoCustomise.m: This sets all parameters individually. We also compare multiple
  instantiations of Add-GP-UCB with different group sizes along with GP-Expected
  Improvement and  DiRect (Dividing Rectangles).

### Citation
If you use this library in your academic work please cite our ICML paper:
"High Dimensional Bayesian Optimisation and Bandits using Additive Models",
Kirthevasan Kandasamy, Jeff Schneider, Barnabas Poczos. International Conference on
Machine Learning, 2015. \\
url: http://www.cs.cmu.edu/~kkandasa/pubs/kandasamyICML15addGPUCB.pdf \\
We use DiRect to optimise the acquisition function. The implementation was taken from
Dan Finkel (2004).

### License
This software is released under GNU GPL v3(>=) License. Please read LICENSE.txt for
more information. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

"Copyright 2015 Kirthevasan Kandasamy"



- For questions/ bug reports please email kandasamy@cs.cmu.edu
