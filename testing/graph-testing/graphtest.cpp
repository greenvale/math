#include <iostream>

#include <graph.hpp>

int main() 
{
    mathlib::Graph<double, double> myGraph;
    
    myGraph.addNode(0, 2.3);
    myGraph.addNode(1, 2.3);
    myGraph.addNode(2, 3.5);
    
    myGraph.addEdge(0, 0, 1, 1.5);
    
    std::cout << myGraph.getEdgeValueById(0) << std::endl;
    
}
