# explainable-opt-code
Code and data for the paper "A Model for Explainable Optimization in Decision Making" by Marc Goerigk and Michael Hartisch, arxiv link to be added.

The code requires Cplex to run.

It provides the following routines to solve shortest path problems:
- Greedy algorithms to find 2 and to find 4 solutions with a decision tree.
- An IP formulation for the same purpose.
- An IP formulation to solve the min-sum-min problem, i.e., to calculate solutions without a decision tree.

To test these methods, a graph generator is provided to create grid graphs. Additionally, real-world data to model Chicago is provided.
