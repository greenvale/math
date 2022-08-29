/*

Graph object class header file
William Denny, 29th August 2022

*/

#pragma once

#include <vector>
#include <assert.h>

namespace mathlib
{
    /* Node object */
    template <class T>
    class Node
    {
    private:
        int m_id;
        T m_value;
    public:
        Node();
        Node(const int& id, const T& value);
        
        int getId();
        T getValue();
        
        void setId(const int& id);
        void setValue(const T& value);       
        
    };
    
    /* Edge object */
    template <class T>
    class Edge 
    {
    private:
        int m_id;
        int* m_nodeIdPair;
        T m_value;
    public:
        Edge();
        Edge(const int& id, const int& nodeId0, const int& nodeId1, const T& value);

        int getId();
        int getNodeId(const int& ref);
        T getValue();
        
        void setId(const int& id);
        void setNodeId(const int& ref, const int& id);
        void setValue(const T& value);   
    };
    
    /* ============================================================== */
    
    /* Node object member functions */ 
    
    template <class T>
    Node<T>::Node() 
    {
    }
    
    template <class T>
    Node<T>::Node(const int& id, const T& value)
    {
        m_id = id;
        m_value = value;
    }
    
    template <class T>
    int Node<T>::getId()
    {
        return m_id;
    }
    
    template <class T>
    T Node<T>::getValue()
    {
        return m_value;
    }
    
    template <class T>
    void Node<T>::setId(const int& id)
    {
        m_id = id;
    }
    
    template <class T>
    void Node<T>::setValue(const T& value)    
    {
        m_value = value;
    }
    
    /* Edge object member functions */
    
    template <class T>
    Edge<T>::Edge() 
    {
    }
    
    template <class T>
    Edge<T>::Edge(const int& id, const int& nodeId0, const int& nodeId1, const T& value)
    {
        m_id = id;
        m_nodeIdPair = new int[2];
        m_nodeIdPair[0] = nodeId0;
        m_nodeIdPair[1] = nodeId1;
        m_value = value;
    }
    
    template <class T>
    int Edge<T>::getId()
    {
        return m_id;
    }
    
    template <class T>
    int Edge<T>::getNodeId(const int& ref)
    {
        return m_nodeIdPair[ref];
    }
    
    template <class T>
    T Edge<T>::getValue()
    {
        return m_value;
    }
    
    template <class T>
    void Edge<T>::setId(const int& id)
    {
        m_id = id;
    }
    
    template <class T>
    void Edge<T>::setNodeId(const int& ref, const int& id)
    {
        m_nodeIdPair[ref] = id;
    }
    
    template <class T>
    void Edge<T>::setValue(const T& value)    
    {
        m_value = value;
    }
    
    
    /* ============================================================== */
    /* ============================================================== */

    /* Graph object */
    template <class T, class U>
    class Graph
    {
    private:
        std::vector<Node<T>> m_nodes;
        std::vector<Edge<U>> m_edges;
    public:
        Graph();
        
        int getNumNodes();
        int getNumEdges();
        std::vector<int> getNodeIds();
        std::vector<int> getEdgeIds();
        
        // ======================================
        // adding and removing nodes
        void addNode(const int& id, const T& value);
        void removeNodeById(const int& id);
        void removeNodesByValue(const T& value);
        
        // getting node indexes
        int getNodeIndexById(const int& id);
        std::vector<int> getNodeIndexesByValue(const T& value);
        
        // getting node(s) id/value by value/id
        std::vector<int> getNodeIdsByValue(const T& value);
        T getNodeValueById(const int& id);
        
        // ======================================
        // adding and removing edges
        void addEdge(const int& id, const int& nodeId1, const int& nodeId2, const U& value);
        void removeEdgeById(const int& id);
        
        // getting edge index by id
        int getEdgeIndexById(const int& id);
        
        // getting edge ids given requirements for edge node pairs
        std::vector<int> getEdgeIdsByNodeIds(std::vector<int> nodeIds);
        
        // get edge value given id
        U getEdgeValueById(const int& id);
    };
    
    /* ============================================================== */

    template <class T, class U>
    Graph<T, U>::Graph()
    {
        m_nodes = {};
        m_edges = {};
    }
    
    template <class T, class U>
    int Graph<T, U>::getNumNodes()
    {
        return m_nodes.size();
    }
    
    template <class T, class U>
    int Graph<T, U>::getNumEdges()
    {
        return m_edges.size();
    }
    
    template <class T, class U>
    std::vector<int> Graph<T, U>::getNodeIds()
    {
        std::vector<int> ids = {};
        int numNodes = getNumNodes();
        for (int i = 0; i < numNodes; ++i)
        {
            ids.push_back(m_nodes[i].getId());
        }
        
        return ids;
    }
    
    template <class T, class U>
    std::vector<int> Graph<T, U>::getEdgeIds()
    {
        std::vector<int> ids = {};
        int numEdges = getNumEdges();
        for (int i = 0; i < numEdges; ++i)
        {
            ids.push_back(m_edges[i].getId());
        }
        
        return ids;
    }
    
    /* ============================================================== */
    
    /* Graph node functionality */
    
    /* adds node with given id and value (errors if id is duplicate of existing id) */
    template <class T, class U>
    void Graph<T, U>::addNode(const int& id, const T& value)
    {
        // check id is not already used in node in m_nodes
        int duplication = 0;
        int numNodes = getNumNodes();
        for (int i = 0; i < numNodes; ++i)
        {
            if (m_nodes[i].getId() == id)
            {
                duplication++;
            }
        }
        assert(duplication == 0);
    
        m_nodes.push_back(Node<T>(id, value));    
    }
    
    template <class T, class U>
    void Graph<T, U>::removeNodeById(const int& id)
    {
        // obtain index and ensure that the node with given id exists in m_nodes
        int index = getNodeIndexById(id);
        assert(index >= 0);
        
        m_nodes.erase(m_nodes.begin() + index);
    }
    
    template <class T, class U>
    void Graph<T, U>::removeNodesByValue(const T& value)
    {
        // obtain vector of indexes of nodes with given value
        std::vector<int> indexes = getNodeIndexesByValue(value);
        
        for (int i = 0; i < indexes.size(); ++i)
        {
            m_nodes.erase(m_nodes.begin() + indexes[i]);
        }
    }
    
    /* returns index of node in m_node with given id (note id duplication is already prevented) */
    template <class T, class U>
    int Graph<T, U>::getNodeIndexById(const int& id)
    {
        int numNodes = getNumNodes();
        for (int i = 0; i < numNodes; ++i)
        {
            if (m_nodes[i].getId() == id)
            {
                return i;
            }
        }
        return -1; // if node with given id not found, return -1
    }
    
    /* returns vector of indexes of nodes in m_nodes with given value */
    template <class T, class U>
    std::vector<int> Graph<T, U>::getNodeIndexesByValue(const T& value)
    {
        std::vector<int> indexes = {};
        
        int numNodes = getNumNodes();
        for (int i = 0; i < numNodes; ++i)
        {
            if (m_nodes[i].getValue() == value)
            {
                indexes.push_back(i);
            }
        }
        
        return indexes;       
    }
    
    /* returns vector of ids for nodes in m_nodes with given value */
    template <class T, class U>
    std::vector<int> Graph<T, U>::getNodeIdsByValue(const T& value)
    {
        std::vector<int> indexes = getNodeIndexesByValue(value);
        std::vector<int> ids = {};
        
        for (int i = 0; i < indexes.size(); ++i)
        {
            ids.push_back(m_nodes[indexes[i]].getId());
        }
        
        return ids;
    }
    
    /* returns value of node with given id (already assumes that there is no id duplication) */
    template <class T, class U>
    T Graph<T, U>::getNodeValueById(const int& id)
    {
        // obtain index and ensure that the node with given id exists in m_nodes
        int index = getNodeIndexById(id);
        assert(index >= 0);
        
        return m_nodes[index].getValue();
    }
    
    /* ============================================================== */
    
    /* Graph edge functionality */
    
    /* add new edge */
    template <class T, class U>
    void Graph<T, U>::addEdge(const int& id, const int& nodeId0, const int& nodeId1, const U& value)
    {
        // check id is not already used in edge in m_edges
        int duplication = 0;
        int numEdges = getNumEdges();
        for (int i = 0; i < numEdges; ++i)
        {
            if (m_edges[i].getId() == id)
            {
                duplication++;
            }
        }
        assert(duplication == 0);
        
        // check node ids provided are valid and correspond to node pair in m_nodes
        int index0 = getNodeIndexById(nodeId0);
        int index1 = getNodeIndexById(nodeId1);
        assert((index0 >= 0) && (index1 >= 0));
        
        m_edges.push_back(Edge<U>(id, nodeId0, nodeId1, value));
    }
    
    template <class T, class U>
    void Graph<T, U>::removeEdgeById(const int& id)
    {
    }
    
    /* returns index of edge in m_edges given id (already assumes no id duplication) */
    template <class T, class U>
    int Graph<T, U>::getEdgeIndexById(const int& id)
    {
        int numEdges = getNumEdges();
        for (int i = 0; i < numEdges; ++i)
        {
            if (m_edges[i].getId() == id)
            {
                return i;
            }
        }
        return -1; // if node with given id not found, return -1
    }
    
    /* returns vector of edge ids that fit node id specs provided (can be a pair or single node id) */
    template <class T, class U>
    std::vector<int> Graph<T, U>::getEdgeIdsByNodeIds(std::vector<int> nodeIds)
    {
        // make sure 1 or 2 ids are given in ids
        assert((nodeIds.size() == 1) || (nodeIds.size() == 2));
        
        std::vector<int> edgeIds = {};
        
        int numEdges = getNumEdges();
        for (int i = 0; i < numEdges; ++i)
        {
            if (nodeIds.size() == 1)
            {
                // if either node ids in edge correspond to given node id, then add edge id to vector
                if ((m_edges[i].getNodeId(0) == nodeIds[0]) || (m_edges[i].getNodeId(1) == nodeIds[0]))
                {
                    edgeIds.push_back(m_edges[i].getId());
                }
            }
            else
            {
                // if node id pair provided, then this pair must be same as that in edge
                if (
                    (m_edges[i].getNodeId(0) == nodeIds[0]) && (m_edges[i].getNodeId(1) == nodeIds[1])
                ||  (m_edges[i].getNodeId(0) == nodeIds[1]) && (m_edges[i].getNodeId(1) == nodeIds[0])
                )
                {
                    edgeIds.push_back(m_edges[i].getId());
                }
            }
        }
        
        return edgeIds;
    }
    
    /* returns value of edge with given id (already assumes that there is no id duplication) */
    template <class T, class U>
    U Graph<T, U>::getEdgeValueById(const int& id)
    {
        // obtain index and ensure that the node with given id exists in m_nodes
        int index = getEdgeIndexById(id);
        assert(index >= 0);
        
        return m_edges[index].getValue();
    }
}
