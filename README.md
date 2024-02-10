# Overview
The repository contains source code for the paper titled [Fast Biconnectivity Restoration in Multi-Robot Systems for Robust Communication Maintenance](https://arxiv.org/pdf/2011.00685.pdf). This paper proposes a graph theoretic algorithm for the _Fast Biconnectivity Restoration_ (FBR) problem. The objective of the FBR problem is to [biconnect](https://en.wikipedia.org/wiki/Biconnected_graph) a barely connected graph (a graph that is connected, but not biconnected) in such a way that the maximum movement of the nodes is minimized. This work utilized the concept of edge augmentation proposed by [Khuller et al.](https://www.sciencedirect.com/science/article/abs/pii/S0196677483710102). To demonstrate the effectiveness of our proposed work, we compare our work with two other existing works on the FBR problem: the [Block Translation Algorithm](https://ieeexplore.ieee.org/document/1316760) proposed by Basu et al. and the [Cascaded Relocation Algorithm](https://ieeexplore.ieee.org/document/4689552) proposed by Abbasi et al. Our proposed algorithm provides a nice tradeoff between accuracy and running time, and outperforms both prior works mentioned above. This repository contains my implementation of the above three papers by Khuller et al., Basu et al., and Abbasi et al.

In this paper, I also develop a QCP (Quadratically Constrained Program) for the FBR problem which can be used to optimally solve the problem. The code to formulate the QCP and solve it is available in this repository. The implementation of the QCP formulation and solution is done using [Gurobi](https://www.gurobi.com/solutions/gurobi-optimizer), a state-of-the-art optimizer for linear and quatratic programs. Although the QCP-based approach solves the FBR problem optimally, it can be used solve only small instances of the FBR problem due to high computational overhead.

In this work, I also test the performance of the proposed algorithm in simulated and real-life environment using a case study. In the case study, a team of robot maintain a biconnected communication topology among them while performing the [persistent monitoring task](https://ieeexplore.ieee.org/abstract/document/8815211). We randomly drop robots to disrupt biconnectivity and make the communication topology barely connected, and then use our proposed solution to biconnect the graph. This case study is performed both in a simulated environment rendered using [OpenGL](https://open.gl/), and also in real-life setting using 6 [CrazyFlie drones](https://www.bitcraze.io/products/crazyflie-2-1/).

# Directory Tree and Code Explanation

* `dataset`: `ds_gen.cpp` contains the code to generate synthetic data for the FBR problem. To this end, I generate a set random points, consider the points as nodes of a graph, and connect the edges between pair of nodes which are at a distance of at most a predefined threshold value, which we call the communication radius. The idea is similar to [Unit Disk Graphs](https://en.wikipedia.org/wiki/Unit_disk_graph). If the communication graph is connected and not biconnected, I add it to the dataset. 
* `algorithms`
  * `ea_scr.cpp`: Our proposed EA-SCR algorithm solves the FBR problem in two phases: Edge Augmentation (EA) and Sequential Cascaded Relocation (SCR). The `ea_scr.cpp` file contains the implementation of this algorithm. The two phases are briefly explained below:
    * In the EA phase, we select a set of edges not present in the graph, which if established, will biconnect the graph. The idea behind the edge augmentation algorithm is similar to the one proposed by [Khuller et al.](https://www.sciencedirect.com/science/article/abs/pii/S0196677483710102). The EA algorithm is based on the concepts of [Block-Cut Tree](https://en.wikipedia.org/wiki/Biconnected_component) and [Minimum Bottleneck Spanning Arborescence](https://en.wikipedia.org/wiki/Minimum_bottleneck_spanning_tree) (MBSA).
    * In the SCR phase, we determine how to move the robots to connect the edges returned by the EA algorithm, such that the maximum movement of a node is minimized. The SCR algorithm is based on the idea of cascaded relocation proposed by [Abbasi et al.](https://ieeexplore.ieee.org/document/4689552). Here, the node pairs indicated by the edges are relocated to make the pairs close enough to reconnect the edges, which might disconnect some of the existing edges. This leads to a cascade of relocations until the whole graph is biconnected.
  * `block_translation.cpp`: This file contains my implementation of the Block Translation Algorithm for the FBR problem proposed by [Basu et al.](https://ieeexplore.ieee.org/document/1316760). The idea is to construct the Block-Cut Tree of the communication graph, and move the blocks as a whole to biconnect the topology. 
  * `cascaded_relocation.cpp`: This file contains my implementation of the Cascaded Relocation Algorithm for the FBR problem proposed by [Abbasi et al.](https://ieeexplore.ieee.org/document/4689552). The idea is to move the nodes minimally to reconnect edges in a recursive manner.
* `qcp-optimizer`
  * `qcp_fbr.cpp`: This code contains the QCP formulation of the FBR problem and solves it using the Gurobi optimizer. The QCP formulation can be found in Section III of the [paper](https://arxiv.org/pdf/2011.00685.pdf). 
  * `qcp_cr.cpp`: This code contains the QCP formulation of the SCR problem and solves it using the Gurobi optimizer. The QCP formulation can be found in Section VI.B of the [paper](https://arxiv.org/pdf/2011.00685.pdf).
* `case-study-pm`: `pm.cpp` contains the code to test the EA-SCR algorithm in the settings of the multi-robot persistent monitoring task. In this task, a team of robots monitor a cluttered environment such that the unobstructed space in the environment is periodically visited by some robot. To this end, the environment is uniformly discretized into cells. The objective of the task is to minimize the total latency of all the unoccupied cells, where latency refers to the time between consequtive visits of a call.
* `visualizer`: `vis.cpp` contains the code to visualize the simulation of the persistent monitoring case study. The visualization is rendered using OpenGL.


# How to Run

1. Generating Dataset: `dataset/ds_gen.cpp` can be used to generate synthetic data for the FBR problem. The data points are generated in a 500 X 500 environment. The number of nodes/robots can be varied using the `dsize` variable on Line 24, and the communication radius can be changed using the `cr` variable on Line 23. The dataset is output in a text file named `data.txt`.
2. Executing Algorithms: `algorithms/ea_scr.cpp` can be used to run the EA-SCR algorithm. The text file containing the dataset must be present in the same directory, and the name of the data file should be updated in the `input_data()` function in the code. The block translation and the cascaded relocation algorithms can be executed in a similar fashion.
3. Finding Optimal Solution: `qcp-optimizer/qcp_fbr.cpp` can be used to find the optimal solution of the FBR problem. The text file containing the dataset must be present in the same directory, and Line 50 in `qcp_fbr.cpp` should be updated to reflect the name of the data file. This code uses the Gorobi optimizer library to perform the optimization. The Gurobi C++ library can be installed with the help of [this tutorial](https://support.gurobi.com/hc/en-us/community/posts/360050201071-How-to-install-gurobi-in-C-Interface).
4. Visualizing Results: `visualizer/vis.cpp` can be used to visualize the persistent monitoring case study. The environment to be visualized can be modified using the variables in Line 22-32 of the code. Here `numobs`, `crange`, `numtargets`, and `lmax` denotes the number of obstacles, communication range, number of targets, and maximum latency, respectively. This code uses the OpenGL API to visualize the simulation. OpenGL can be installed with the help of [this tutorial](https://www.opengl-tutorial.org/beginners-tutorials/tutorial-1-opening-a-window/).


# Links

## Papers
* [Fast Biconnectivity Restoration in Multi-Robot Systems for Robust Communication Maintenance](https://arxiv.org/pdf/2011.00685.pdf)
* [Approximation Algorithms for Graph Augmentation](https://www.sciencedirect.com/science/article/abs/pii/S0196677483710102)
* [Movement control algorithms for realization of fault-tolerant ad hoc robot networks](https://ieeexplore.ieee.org/document/1316760)
* [Movement-Assisted Connectivity Restoration in Wireless Sensor and Actor Networks](https://ieeexplore.ieee.org/document/4689552)

## Tools
* [Gurobi](https://www.gurobi.com/solutions/gurobi-optimizer)
* [OpenGL](https://open.gl/)
* [CrazyFlie](https://www.bitcraze.io/products/crazyflie-2-1/)

