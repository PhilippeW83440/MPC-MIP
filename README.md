

# AA222_Project
Final project on Trajectory Optimization

## Elastic Model for Collision Avoidance  
Sometimes a constraint can not be satisfied, but we would like nevertheless to return a solution that minimizes deviation to our constraint. Instead of just throwing INFEASIBLE ...  
For example, if we can not maintain a safety distance of 4 meters, then it may be OK to use 3 meters ..  
=> Allow to violate constraints at a cost. Introduce slack variables and add them to the objective as penalty.  

![\begin{split}\begin{array}{lr}   \minimize    & c^Tx, \\   \st          & a^Tx\leq b. \end{array}\end{split}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Bsplit%7D%5Cbegin%7Barray%7D%7Blr%7D%20%20%20%5Cminimize%20%20%20%20%26%20c%5ETx%2C%20%5C%5C%20%20%20%5Cst%20%20%20%20%20%20%20%20%20%20%26%20a%5ETx%5Cleq%20b.%20%5Cend%7Barray%7D%5Cend%7Bsplit%7D)

which might be causing infeasibility. Then create a new variable y and form the problem which contains:

![\begin{split}\begin{array}{lr} \minimize    & c^Tx+y, \\ \st          & a^Tx\leq b+y. \end{array}\end{split}](https://render.githubusercontent.com/render/math?math=%5Cbegin%7Bsplit%7D%5Cbegin%7Barray%7D%7Blr%7D%20%5Cminimize%20%20%20%20%26%20c%5ETx%2By%2C%20%5C%5C%20%5Cst%20%20%20%20%20%20%20%20%20%20%26%20a%5ETx%5Cleq%20b%2By.%20%5Cend%7Barray%7D%5Cend%7Bsplit%7D)  
Solving this problem will reveal by how much the constraint needs to be relaxed in order to become feasible.   
This is equivalent to inspecting the infeasibility certificate but may be more intuitive.

## Optimization Algorithms  
* [Log Barrier with equality constraints](https://www.stat.cmu.edu/~ryantibs/convexopt/lectures/barr-method.pdf)  
* [Equality constrained Newton: Stephen Boyd book p.526,545/546,548,552/553](https://web.stanford.edu/~boyd/cvxbook/bv_cvxbook.pdf)    
* [Primal-Dual Interior-Point Methods](https://www.stat.cmu.edu/~ryantibs/convexopt/lectures/primal-dual.pdf)  
* [Augmented Lagrangian](http://www.cs.cmu.edu/~pradeepr/convexopt/Lecture_Slides/Augmented-lagrangian.pdf) and class book section 10.8 p.183  
* [Disjunctive Inequalities](https://optimization.mccormick.northwestern.edu/index.php/Disjunctive_inequalities)  
* [Ipopt paper](http://cepac.cheme.cmu.edu/pasilectures/biegler/ipopt.pdf)  
  - Equality Constraints in Ipopt: Cf section 3.5.  Handling Problems Without a Strict Relative Interior


## Trajectory Optimization Resources  
* [Trajectory Optimization](http://www.matthewpeterkelly.com/tutorials/trajectoryOptimization/index.html)  
* [Handling disjunctive behaviours: Practical Methods for Optimal Control and Estimation Using NLP, Ch 1.15](https://books.google.fr/books?id=n9hLriD8Lb8C&printsec=frontcover#v=onepage&q&f=false)  
* [Speed Profile Planning in Dynamic Environments via Temporal Optimization](http://www.cs.cmu.edu/~cliu6/files/iv17-2.pdf)
* [Speed Profile Planning Poster](http://www.cs.cmu.edu/~cliu6/files/iv17-2poster.pdf), [Dissertation ch 5.4 & 6](http://www.cs.cmu.edu/~cliu6/files/dissertation.pdf) and [code](https://github.com/changliuliu/TemporalOptimization)   
* [The Convex Feasible Set algorithm for real time optimization in Motion Planning](https://arxiv.org/pdf/1709.00627.pdf) 
* [Autonomous Vehicle Control: A Nonconvex Approach for Obstacle Avoidance](https://www.researchgate.net/publication/303905926_Autonomous_Vehicle_Control_A_Nonconvex_Approach_for_Obstacle_Avoidance)  
* [Explicit MPC (slides 15-16)](https://stanford.edu/class/ee364b/lectures/mpc_slides.pdf) and computation of Invariant Sets with [MPT3](https://www.mpt3.org/UI/Invariance)  
* [Fast MPC (Explicit MPC)](https://web.stanford.edu/~boyd/papers/pdf/fast_mpc.pdf)  
* [Optimization-Based Collision Avoidance](https://arxiv.org/pdf/1711.03449.pdf)




## Path2d generation  
* [Generating paths with polynomial spirals](https://github.com/qiaoxu123/Self-Driving-Cars/blob/master/Part4-Motion_Planning_for_Self-Driving_Cars/Module7-Putting_it_all_together-Smooth_Local_Planning/Module7-Putting_it_all_together-Smooth_Local_Planning.md), [CMU DARPA paper](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.225.4980&rep=rep1&type=pdf)  
* [Path model](https://www.ri.cmu.edu/pub_files/2012/5/ICRA12_xuwd_Final.pdf)
