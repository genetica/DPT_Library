/******************************************************************************
Discrete Pulse Transform Library.
DPT Library v0.1

Copyright (C) 2013, Gene Stoltz
All rights reserved.
******************************************************************************/

#include <iostream>
#include <fstream>

#include <cstdlib>

#include <time.h>

//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/core/core.hpp>
//#include <opencv2/opencv.hpp>

#include "DPT.h"

using namespace std;

WGraphNode *zero;

    //static clock_t tot = 0;

// include sort then use
//#define USE_SORT

// qsort function for compare of WGraph Node pointers
int compareWGraphNodeP (const void * TNode1, const void * TNode2)
{
    //cout << "\nT1 " << *(WGraphNode **)TNode1 << endl;
    //cout << "\nT2 " << *(WGraphNode **)TNode2 << endl;

    return (*(WGraphNode **)TNode1 - *(WGraphNode **)TNode2);
}
// Remove double connections and sort connections
inline void DeleteDoubleConnections(WGraphNode *&CNode)
{
    // sorteer array
    qsort(CNode->ptr_connected,CNode->n_connected,sizeof(CNode->ptr_connected),compareWGraphNodeP);

    // remove doubles
    int index = 0;
    for (int i = 1;  i < CNode->n_connected; ++i)
    {
        if (CNode->ptr_connected[i] != CNode->ptr_connected[index])
        {
            //indien verskil, las element i by
            ++index;
            CNode->ptr_connected[index] = CNode->ptr_connected[i];
        }
    }
    CNode->n_connected = index + 1;
}

// C++ I/O
int CreateGraph(int n_nodes, int *&connectivity, int *&data_structure, WGraphNode **&index, WGraphNode *&node, PGraphNode *&base)
{
    // Add extra point for boundary conditions
    ++n_nodes;

    // create index
    index = new WGraphNode * [n_nodes];
    // create node memory
    node = new WGraphNode [n_nodes];
    // create pulse base - maximum number of pulses can not exceed the amount of starting
    //                     nodes and the pulse graph contains the datapoint positions
    base = new PGraphNode [2*n_nodes];

    node[0].empty = 0;
    node[0].ptr_merged = 0;
    node[0].size = n_nodes + 1;
    node[0].value = 0;
    node[0].n_construct = 0;
    node[0].ptr_connected = new WGraphNode *[1];
    node[0].ptr_construct = new PGraphNode *[1];
    node[0].n_connected = 0;
    node[0].index = 0;
    index[0] = (&node[0]);
    #ifdef DATA_RANGE
    node[0].histogram = new int [DATA_RANGE];
    #endif
    zero = index[0];


        int n_dim = connectivity[0];
        // position of #connections - see connectivity data structure
        int p_conn = n_dim + 1;
        // number of connections
        int n_conn = connectivity[p_conn];
        int *coordinates = new int[n_dim];
        // addition is the additional scaling required to increase one point on axis
        int *addition = new int[n_dim];
        addition[0] = 1;
        for (int i = 1; i < n_dim; ++i)
        {
            addition[i] = addition[i - 1] * connectivity[1 + (i - 1)];
        }
 //   int tel = 1; // start tel at 1 to skip 0-inf node
//    while (!data.eof())
    for (int tel = 1; tel < n_nodes; ++tel)
    {
//        data >> value;

        // save memory adres in index
        index[tel] = (&node[tel]);

        // initialize nodes in work graph
        node[tel].empty = 0;
        node[tel].size = 1;
        node[tel].value = data_structure[tel - 1];   // value
        node[tel].n_construct = 1;
        node[tel].ptr_construct = new PGraphNode * [1];
        node[tel].ptr_construct[0] = &base[tel - 1];
        node[tel].index = tel;
        node[tel].ptr_merged = 0;
        //create histogram for working graph
        #ifdef DATA_RANGE

            node[tel].histogram = new int [DATA_RANGE];

            for (int k = 0; k < DATA_RANGE; k++)
                node[tel].histogram[k] = 0;
            node[tel].histogram[node[tel].value] = 1;
        #endif

        // initialize base nodes in pulse graph
        // use tel - 1 as the pulse graph does not contain 0-inf node
        base[tel - 1].size = 1;
        base[tel - 1].height = 0;
        base[tel - 1].n_construct = tel - 1; // position in datastructure

        base[tel - 1].Connected_Pulse = 0;
        //base[tel - 1].status = 0;
        //base[tel - 1].CLand = 0;
        base[tel - 1].index = tel - 1;

        // create histotram for base nodes
        #ifdef DATA_RANGE

            base[tel - 1].histogram = new int [DATA_RANGE];

            for (int k = 0; k < DATA_RANGE; ++k)
                base[tel - 1].histogram[k] = 0;
            base[tel - 1].histogram[node[tel].value] = 1;
        #endif


        // graph connectivity - define by using offset
        node[tel].n_connected = n_conn; // connectivity saved in pos of 0 element
        node[tel].ptr_connected = new WGraphNode * [node[tel].n_connected];
        node[tel].ptr_connected[0] = 0;

/*  Connectivity Data structure
        #Dimensions - n
        Dimension[1]
        ...
        Dimension[n]
        #Connections - k     (p_conn)
        Delta[0][0]
        ...
        Delta[0][n]
        ...
        Delta[k][0]
        ...
        Delta[k][n]
*/
        // find co-ordinates in data space ([0..k])
        int place = tel - 1;
        for (int dim = n_dim - 1; dim >= 0; --dim)
        {
    //        if (place > addition[dim])
     //       {
                int mod_pos = place%addition[dim];
                coordinates[dim] = (place - mod_pos)/addition[dim];
                place = mod_pos;
      //      }
      //      else
            {
       //         coordinates[dim] = 0;
            }
        }
  //      if (tel < addition[1])
  //          coordinates[0] = tel - 1;
  //      cout << "coord" << coordinates[0] << " " << coordinates[1] << endl;

        int c_conn = p_conn;
        // set connections
//cout << tel << "\t";
        for (int j = 0; j < node[tel].n_connected ; ++j)
        {
//cout << "\nconnection " << j << "\n";
            int new_pos;
            place = tel;
            for (int dim = 0; dim < n_dim; ++dim)
            {
                ++c_conn;
                new_pos = coordinates[dim] + connectivity[c_conn];
//cout << " pos " << coordinates[dim] << " new " << new_pos << " con " << connectivity[c_conn] << "\n";
                if ((new_pos < 0)||(new_pos >= connectivity[dim + 1]))
                {
                    place = 0;
                    // stop for - connection already out of bounds
                    c_conn += n_dim-dim-1;
                    dim = n_dim;
                }
                else
                {
                    place += connectivity[c_conn]*addition[dim];
                }
            }
//cout << " p " << place << "\t";
            node[tel].ptr_connected[j] = &node[place];
        }
//        cout << endl;
/*
        WGraphNode *CNode = &(node[tel]);
        qsort(CNode->ptr_connected,CNode->n_connected,sizeof(CNode->ptr_connected),compareWGraphNodeP);

        cout << tel << "\t" << CNode << "\t";
        for (int j = 0; j < node[tel].n_connected ; ++j)
        {
            cout << node[tel].ptr_connected[j] << " ";
        }
        cout << endl;
        getchar();
*/
    }

/*
    for (int tel = 0; tel < n_nodes; ++tel)
    {
        cout << tel << "\t"  << index[tel] << " ";
        for (int j = 0; j < index[tel]->n_connected ; ++j)
        {
            cout << index[tel]->ptr_connected[j] << " ";
        }
        cout << index[tel]->value << endl;
    }
    getchar();
*/
    return n_nodes;
}

#ifdef DATA_RANGE
int UpdateHistogram(WGraphNode *&CNode, WGraphNode *&TNode)
{
    for (int i = 0; i < DATA_RANGE; i++)
    {
        TNode->histogram[i] += CNode->histogram[i];
    }
    delete [] CNode->histogram;

    return 0;
}
#endif

inline int UpdateConnectedComp(WGraphNode *&CNode, WGraphNode *&TNode)
{

    //clock_t begin;
    //clock_t end;
    //begin = clock();

    //int k = 0;
    int new_n_connected;

    // Update connected Components in Test Node
    WGraphNode **temp_WGraphNode;


    // create temp variable
    new_n_connected = TNode->n_connected + CNode->n_connected; // mius one, remove zero connection
    temp_WGraphNode = new WGraphNode * [new_n_connected];

    int counter_TNode = 0;
    int counter_CNode = 0;
    int counter = 0;
    /*
cout << "start\n" << "TNode " << TNode << endl;
    for (int n = 0; n < TNode->n_connected; ++n)
    {
        cout << TNode->ptr_connected[n] << "\t";
    }
cout << endl << "CNode " << CNode << endl;
    for (int n = 0; n < CNode->n_connected; ++n)
    {
        cout << CNode->ptr_connected[n] << "\t";
    }
cout << endl << "after\n";
*/
    while ((counter_TNode < TNode->n_connected)||(counter_CNode < CNode->n_connected))
    {
        // If connection in TNode equal to CNode skip connection
        if (TNode->ptr_connected[counter_TNode] == CNode)
            ++counter_TNode;
        // If connection in CNode equal to TNode skip connection
        if (CNode->ptr_connected[counter_CNode] == TNode)
            ++counter_CNode;

        // if both nodes is still within assigned memory
        if ((counter_TNode < TNode->n_connected)&&(counter_CNode < CNode->n_connected))
        {
            // if CNode connection is smaller than TNode connection - add CNode connection
            if (CNode->ptr_connected[counter_CNode] < TNode->ptr_connected[counter_TNode])
            {
                temp_WGraphNode[counter] = CNode->ptr_connected[counter_CNode];
                ++counter;
                ++counter_CNode;
            }
            else
            // if TNode connection is smaller than CNode connection - add TNode connection
            if (CNode->ptr_connected[counter_CNode] > TNode->ptr_connected[counter_TNode])
            {
                temp_WGraphNode[counter] = TNode->ptr_connected[counter_TNode];
                ++counter;
                ++counter_TNode;
            }
            else
            // if TNode  connection and CNode connection is equal - add only CNode
            // if (CNode->ptr_connected[counter_CNode] == TNode->ptr_connected[counter_TNode])
            {
                if (counter == 0)
                {
                    temp_WGraphNode[counter] = CNode->ptr_connected[counter_CNode];
                    ++counter;
                }
                else
                if (temp_WGraphNode[counter - 1] != CNode->ptr_connected[counter_CNode])
                {
                    temp_WGraphNode[counter] = CNode->ptr_connected[counter_CNode];
                    ++counter;
                }

                ++counter_CNode;
                ++counter_TNode;
            }
        }
        else
        // if TNode counter is out of bounds - add all CNodes connections
        if (counter_TNode < TNode->n_connected)
        {
            temp_WGraphNode[counter] = TNode->ptr_connected[counter_TNode];
            ++counter;
            ++counter_TNode;
        }
        else
        // if CNode counter is out of bounds - add all TNode connections
        if (counter_CNode < CNode->n_connected)
        {
            temp_WGraphNode[counter] = CNode->ptr_connected[counter_CNode];
            ++counter;
            ++counter_CNode;
        }
    }
    // if no connections exist assign all other Nodes.
    if (TNode->n_connected == 0)
    {
        for (int n = 0; n < CNode->n_connected; ++n)
            temp_WGraphNode[n] = CNode->ptr_connected[n];
        counter = CNode->n_connected;
    }
    else
    if (CNode->n_connected == 0)
    {
        for (int n = 0; n < TNode->n_connected; ++n)
            temp_WGraphNode[n] = TNode->ptr_connected[n];
        counter = TNode->n_connected;
    }
/*
    // copy test node connected components to temp
    for (int n = 0; n < TNode->n_connected; ++n)
    {
        if ((TNode->ptr_connected[n] == CNode)||(TNode->ptr_connected[n] == TNode))
        {
            --k;
        }
        else
        {
            temp_WGraphNode[k] = TNode->ptr_connected[n];
        }
        ++k;
    }

    // copy current nodes connected components to temp
    for (int n = 0; n < CNode->n_connected; ++n)
    {
        if ((CNode->ptr_connected[n] == TNode)||(CNode->ptr_connected[n] == CNode))
        {
            --k;
        }
        else
        {
            temp_WGraphNode[k] = CNode->ptr_connected[n];
        }
        ++k;
    }
*/
    // update new node connected component information
    TNode->n_connected = counter;

    // delete old connected comp
    delete [] TNode->ptr_connected;

    // Assign new connected comp pointer
    TNode->ptr_connected = temp_WGraphNode;

        /* // No need to copy just assign pointer
        TNode->ptr_connected = new WGraphNode * [TNode->n_connected];
        // copy all components to new node
        for (int n = 0; n < TNode->n_connected; n++)
        {
            TNode->ptr_connected[n] = temp_WGraphNode[n];
        }
        // delete temporary data
        delete [] temp_WGraphNode;
        */
/*
    for (int n = 0; n < TNode->n_connected; ++n)
    {
        cout << TNode->ptr_connected[n] << "\t";
    }
cout << endl;
getchar();
*/
    // Delete merged node connected comp
    delete [] CNode->ptr_connected;
    CNode->ptr_connected = 0;
    // increase size of new node
    TNode->size += CNode->size;

#if defined USE_SORT
    DeleteDoubleConnections(TNode);
#endif

    //end = clock();
    //tot += (end - begin);
    return 0;
}

inline int UpdateConstructionComp(WGraphNode *&CNode, WGraphNode *&TNode)
{
    int new_n_construct;

    PGraphNode **temp_PGraphNode;

    new_n_construct = TNode->n_construct + CNode->n_construct;

    temp_PGraphNode = new PGraphNode * [new_n_construct];

    // copy construct components to temp variable
    for (int n = 0; n < TNode->n_construct; ++n)
    {
        temp_PGraphNode[n] = TNode->ptr_construct[n];
    }

    delete [] TNode->ptr_construct;
    // change pointer array
    TNode->ptr_construct = temp_PGraphNode;

        /*
        TNode->ptr_construct = new PGraphNode * [new_n_construct];
        // copy construction components from temp variable
        for (int n = 0; n < TNode->n_construct; n++)
        {
            TNode->ptr_construct[n] = temp_PGraphNode[n];
        }

        delete [] temp_PGraphNode;\
        */

    // add current node construction components to test node
    for (int n = 0; n < CNode->n_construct; ++n)
    {
        TNode->ptr_construct[TNode->n_construct + n] = CNode->ptr_construct[n];
    }

    TNode->n_construct = new_n_construct;

    delete [] CNode->ptr_construct;

    return 0;
}

inline bool CheckNodeConnectedness(WGraphNode  *&CNode, int &i_component)
{
    bool removed = false;
    // check if connected node is empty - update until connected is not empty
    WGraphNode *TNode = CNode->ptr_connected[i_component];

    while (TNode->empty == 1)
    {
        // get empty node merged adres
        TNode = TNode->ptr_merged;
        CNode->ptr_connected[i_component] = TNode;
        // check if connected node references itself
    }

    // if TNode and CNOde is equal remove TNode connection
    if (CNode == TNode)
    {
        --CNode->n_connected;
        // if node references itself - delete connected node from connected components
        for (int k = i_component; k < CNode->n_connected; ++k)
        {
            CNode->ptr_connected[k] = CNode->ptr_connected[k+1];
        }
        // check if i_component is out off bounds
        removed = true;
    }

    return removed;
}

inline int CheckFeature(int &n_nodes, int &feature_count,WGraphNode *&CNode, WGraphNode **&feature)
{
    // output 1 - if feature otherwise 0
    // also update current feature index

    int feature_test = 0;
    int feature_return = 0;
    int start_index = (CNode->ptr_connected[0] == 0)?1:0;

    // stop boundary condition of becoming a feature
    if (CNode->size <= n_nodes)
    {
        // check if node is max/min feature or nothing
        for (int m = start_index; m < CNode->n_connected; ++m)
        {
            if (CNode->value < CNode->ptr_connected[m]->value)
                --feature_test;
            else
                ++feature_test;
        }
        // check result
        if (feature_test == -1*CNode->n_connected + start_index)
        {
            // min feature - add to feature list
            feature[feature_count] = CNode;
            CNode->empty = 2;       // show feature is min
            feature_return = 1;
        }
        else
        if (feature_test == CNode->n_connected - start_index)
        {
            //max feature -- add to feature list
            feature[feature_count] = CNode;
            CNode->empty = 3;       // 3 - show feature is max
            feature_return = 1;
        }
    }
    return feature_return;
}

// used to minimize nodes on graph before DPT and extract current features
 int OptimizeGraph(int n_nodes, WGraphNode **&feature, WGraphNode **&index, WGraphNode *&node)
{
    // returns the number of features
    int number_of_features = 0;

    feature = new WGraphNode * [n_nodes];

    int node_count = n_nodes;

    WGraphNode  *CNode,  // Current Node
                *TNode;  // Test Node

    for (int j = 1; j < node_count; ++j)
    {
//        cout << j << "\t" << node_count << endl;
        CNode = index[j];

        // check all connections and compare for possible merger
        // check if zero node exist in node, this determines starting index
        for (int m = 0; m < CNode->n_connected; ++m)
        {
            // check if connected node has merged with another - if yes update connection
            // if node is removed recheck node
            if (CheckNodeConnectedness(CNode,m))
            {
                --m;
            }
            else
            {
                TNode = CNode->ptr_connected[m];

                // check connected comp same values
                if (CNode->value == TNode->value)
                {
                    // update connected components of TNode
                    for (int t = 0;t < TNode->n_connected; ++t)
                    {
                        CheckNodeConnectedness(TNode,t);
                    }

            // Update connected Components in Test Node
                    UpdateConnectedComp(CNode, TNode);

            // Update construction Components in Test Node
                    UpdateConstructionComp(CNode,TNode);

            // Update Histogram in Test Node
                    #ifdef DATA_RANGE
                        UpdateHistogram(CNode,TNode);
                    #endif

                    // set current node to empty
                    CNode->empty = 1;
                    CNode->ptr_merged = TNode;
                    // stop merging this pulse
                    m = CNode->n_connected;

                    // remove node from index
                    index[CNode->index] = index[node_count-1];  //       index[j] = index[i-1];
                    index[node_count-1]->index = CNode->index;  //  index[i-1]->index = j;

                    index[node_count-1] = CNode;  //       index[j] = index[i-1];
                    index[node_count-1]->index = node_count-1;  //  index[i-1]->index = j;

                    // reduce total number of nodes in index
                    --node_count;
                    // redo current j index
                    --j;
                }
            }
        }
    }


/*
    for (int j = 0; j < n_nodes; ++j)
    {
        cout << j << "\t" << &(node[j]) << " " << node[j].ptr_merged << " " << node[j].empty;
        if (node[j].ptr_merged == 0)
        {
            cout << "\n\t" << node[j].n_connected << " " << node[j].value << " " << node[j].size << endl;
            cout << " CONNECTED \n";
        for (int m = 0; m < node[j].n_connected; ++m)
        {
            cout << node[j].ptr_connected[m] << " " << node[j].empty << " ";
        }
        cout << endl;

        }
        cout << endl;
    }
    getchar();
*/
    // update and smooth all nodes not merged
    for (int j = 1; j < n_nodes; ++j)
    {
        if (node[j].ptr_merged == 0)
        {
            CNode = &(node[j]);
            for (int m = 0; m < node[j].n_connected; ++m)
            {
                if (CheckNodeConnectedness(CNode,m))
                    --m;
            }
            DeleteDoubleConnections(CNode);
        }
    }
    // extract features
    for (int j = 1; j < n_nodes; ++j)
    {
//        cout << j << "\t" << n_nodes << endl;
        CNode = &(node[j]);
        if (CNode->ptr_merged == 0)
        {
            number_of_features += CheckFeature(n_nodes,number_of_features,CNode,feature);
        }
    }

    //cout << "number of features " << number_of_features << endl;

    return number_of_features;
}

//called within DPT
inline int OptimizeWorkGraph(WGraphNode *&CNode)
{
    WGraphNode *TNode;
    int stop = 0;
    int tel  = 0;

    // run thourgh connected components
    while ((stop ==0)&&(CNode->n_connected > 0))
    {

        CheckNodeConnectedness(CNode,tel);

        TNode = CNode->ptr_connected[tel];

        // check connected comp same values
        if (CNode->value == TNode->value)
        {
            // Update connected Components in Test Node
            UpdateConnectedComp(CNode, TNode);

            // Update construction Components in Test Node
            UpdateConstructionComp(CNode,TNode);

            // Update Histogram in Test Node
            #ifdef DATA_RANGE
            UpdateHistogram(CNode,TNode);
            #endif

            // set current node to empty
            CNode->empty = 1;
            CNode->ptr_merged = TNode;

            // Check TNode connected components now
            CNode = TNode;
            tel = -1;
        }

        ++tel;
        // stop if all connected components was tested for equal values
        // and none matched
        if (tel >= CNode->n_connected)
        {
            stop = 1;
        }
    }
    return 0;
}

int AddPulse(int output_value, int &pulse_index,WGraphNode *&CNode,WGraphNode *&TNode,PGraphNode *&pulse )
{

    ++pulse_index;
    pulse[pulse_index].height = CNode->value - output_value;
    pulse[pulse_index].size   = CNode->size;
    pulse[pulse_index].Connected_Pulse = 0;
    pulse[pulse_index].n_construct = CNode->n_construct;
    pulse[pulse_index].ptr_construct = new PGraphNode * [pulse[pulse_index].n_construct];
	//pulse[pulse_index].status = 0;
	//pulse[pulse_index].CLand = 0;
	pulse[pulse_index].index = pulse_index;

    //copy construction
    for (int j = 0; j < pulse[pulse_index].n_construct; j++)
    {
        // copy node construction to pulse
        pulse[pulse_index].ptr_construct[j] = CNode->ptr_construct[j];
        // update current pulse construction pulses
        pulse[pulse_index].ptr_construct[j]->Connected_Pulse = &pulse[pulse_index];
    }

    // Copy Histogram in Test Node
    #ifdef DATA_RANGE
        pulse[pulse_index].histogram = new int [DATA_RANGE];

        for (int i = 0; i < DATA_RANGE; i++)
        {
            pulse[pulse_index].histogram[i] = CNode->histogram[i];
        }
    #endif

    //Remove current node by combining to TNode
    //Update Connected Components
    UpdateConnectedComp(CNode, TNode);              // TNode size update here as well.

    // Update Construction Components
    CNode->n_construct = 1;
    delete [] CNode->ptr_construct;                 // Remove all construction nodes as evrything
    CNode->ptr_construct = new PGraphNode * [1];    // is summed in the pulse created.
    CNode->ptr_construct[0] = &pulse[pulse_index];  // Add only pulse to construction

    UpdateConstructionComp(CNode,TNode);

    // Update Histogram in Test Node
    #ifdef DATA_RANGE
        UpdateHistogram(CNode,TNode);
    #endif

    // set current node to empty
    CNode->empty = 1;
    CNode->ptr_merged = TNode;

//    cout << "Pulse :" << CNode->value << "\t"<< pulse[pulse_index].height << "\t" << pulse[pulse_index].size << endl;
    return 0;
}

int DPT_decomposition(int &number_of_features, int &n_nodes, WGraphNode **&feature, PGraphNode *&pulse, WGraphNode **&index)
{
    int output_value = 0;

    int pulse_index = n_nodes-1;

    WGraphNode *CNode,
               *TNode,
               *tempNode;

    int n = 0;
    int min = 1;

    int telhoev = 0;
    int tempSize;
    //clock_t begin;
    //clock_t end1;
    //clock_t end2;

    int feature_i = 0;

    bool updated_feature = false;

    while (number_of_features > 0)
    {
        n = min;
//        begin = clock();


/*
cout << endl << n <<endl;
//if (n > 20000)
{
for (int i = 0;i < number_of_features; i++)
{
    int teller = 0;
    for (int k = 0; k < feature[i]->n_connected; k++)
        for(int y = k + 1; y < feature[i]->n_connected; y++)
        {
            if (feature[i]->ptr_connected[y] == feature[i]->ptr_connected[k])
            {
                feature[i]->ptr_connected[y] = feature[i]->ptr_connected[feature[i]->n_connected - 1];
                feature[i]->n_connected--;
                teller++;
            }
        }

//    cout << "bloom " << teller << " ";
    cout << i <<"\t";
    cout << feature[i]->value << "\t";
    cout << feature[i]->size  << "\t";
    cout << feature[i]->empty  << "\t";
    cout << feature[i]->n_connected << "\t";
    cout << feature[i]->n_construct << endl;

}
getchar();
}
*/

       min = 2*n_nodes;
// min max
        feature_i = 0;
        for (int i = 0;i < number_of_features; ++i)
        {
            CNode = feature[i];

            if (CNode->empty == 1)
            {
                --number_of_features;
                // if not a feature - delete from index
                feature[i] = feature[number_of_features];
                --i;
            }
            else
            if (CNode->size <= n)
            {
                updated_feature = false;

                // check all connections and compare for possible merger
                for (int m = 0; m < CNode->n_connected; ++m)
                {
                    // check if connected node has merged with another - if yes update connection
                    // if node is removed recheck node
                    if (CheckNodeConnectedness(CNode,m))
                    {
                        --m;
                    }
                    else
                    {
                        TNode = CNode->ptr_connected[m];

                        // check connected comp same values
                        if (CNode->value == TNode->value)
                        {
                            for (int t = 0;t < TNode->n_connected; ++t)
                            {
                                if (CheckNodeConnectedness(TNode,t))
                                    --t;
                            }

                    // Update connected Components in Test Node
                            UpdateConnectedComp(CNode, TNode);

                    // Update construction Components in Test Node
                            UpdateConstructionComp(CNode,TNode);

                    // Update Histogram in Test Node
                            #ifdef DATA_RANGE
                                UpdateHistogram(CNode,TNode);
                            #endif

                            // set current node to empty
                            CNode->empty = 1;
                            CNode->ptr_merged = TNode;
                            // stop merging this pulse
                            m = CNode->n_connected;

                            // Update possible feature
                            feature[i] = TNode;
                            updated_feature = true;

                            // redo current feature i index
                            --i;
                        }
                    }
                }

                if (updated_feature == false)
                {
                    if (CheckFeature(n_nodes,i,CNode,feature) == 0)
                    {
                        --number_of_features;
                        // if not a feature - delete from index
                        feature[i] = feature[number_of_features];
                        // recheck current index
                        --i;
                    }
                    else
                    {
                        if ((CNode->size == n))
                        {
                            if ( CNode->empty == 2)
                            {


                            // find component with closest value
                            output_value = 1000000;
                            for (int j = 0; j < CNode->n_connected; ++j)
                            {
                                if (CNode->ptr_connected[j]->value < output_value)
                                {
                                    output_value = CNode->ptr_connected[j]->value;
                                    //save node with required value
                                    TNode = CNode->ptr_connected[j];
                                }
                            }
                            // update TNode connected components
                            for (int t = 0; t < TNode->n_connected; ++t)
                            {
                                if (CheckNodeConnectedness(TNode,t))
                                    --t;
                            }

                            //Create Pulse
                            AddPulse(TNode->value,pulse_index,CNode,TNode,pulse);

                            CNode = TNode;
                            feature[i] = CNode;
                            }
                            else
                            // if (CNode->empty == 3)
                            {
                                //move feature to beginning of feature table
                                tempNode = feature[feature_i];
                                feature[feature_i] = feature[i];
                                feature[i] = tempNode;
                                ++feature_i;
                            }
                        }
                    }
                }
            }
            if ((i>=0))
            {
                tempSize = feature[i]->size;
                if ((tempSize < min)&&(tempSize != n))
                    min = tempSize;
            }


        }
//end1 = clock();

// max min
        for (int i = 0;i < feature_i; ++i)
        {
            CNode = feature[i];

            if (CNode->empty == 1)
            {
                --number_of_features;
                // if not a feature - delete from index
                feature[i] = feature[number_of_features];
                --i;
            }
            else
            if (CNode->size <= n)
            {
                updated_feature = false;

                // check all connections and compare for possible merger
                // check if zero node exist in node, this determines starting index
                for (int m = 0; m < CNode->n_connected; ++m)
                {
                    // check if connected node has merged with another - if yes update connection
                    // if node is removed recheck node
                    if (CheckNodeConnectedness(CNode,m))
                    {
                        --m;
                    }
                    else
                    {
                        TNode = CNode->ptr_connected[m];

                        // check connected comp same values
                        if (CNode->value == TNode->value)
                        {
                            for (int t = 0;t < TNode->n_connected; ++t)
                            {
                                if (CheckNodeConnectedness(TNode,t))
                                    --t;
                            }

                    // Update connected Components in Test Node
                            UpdateConnectedComp(CNode, TNode);

                    // Update construction Components in Test Node
                            UpdateConstructionComp(CNode,TNode);

                    // Update Histogram in Test Node
                            #ifdef DATA_RANGE
                                UpdateHistogram(CNode,TNode);
                            #endif

                            // set current node to empty
                            CNode->empty = 1;
                            CNode->ptr_merged = TNode;
                            // stop merging this pulse
                            m = CNode->n_connected;

                            // Update possible feature
                            feature[i] = TNode;
                            updated_feature = true;

                            // redo current feature i index
                            --i;
                        }
                    }
                }

                if (updated_feature == false)
                {
                    if (CheckFeature(n_nodes,i,CNode,feature) == 0)
                    {
                        --number_of_features;
                        // if not a feature - delete from index
                        feature[i] = feature[number_of_features];
                        // recheck current index
                        --i;
                    }
                    else
                    {
                        if ((CNode->empty == 3)&&(CNode->size == n))
                        {
                            // find component with closest value
                            output_value = -1000000;
//cout << endl << CNode << " " << CNode->value << " ";
                            for (int j = 0; j < CNode->n_connected; ++j)
                            {
//cout << CNode->ptr_connected[j] << " " << CNode->ptr_connected[j]->value << " ";
                                if (CNode->ptr_connected[j]->value > output_value)
                                {
                                    output_value = CNode->ptr_connected[j]->value;
                                    //save node with required value
                                    TNode = CNode->ptr_connected[j];
                                }
                            }
                            // update TNode connected components
                            for (int t = 0; t < TNode->n_connected; ++t)
                            {
                                if (CheckNodeConnectedness(TNode,t))
                                    --t;
                            }

                            //Create Pulse
                            AddPulse(TNode->value,pulse_index,CNode,TNode,pulse);

                            CNode = TNode;
                            feature[i] = CNode;
                        }
                    }
                }
            }
            if ((i>=0)&&(feature[i]->size < min))
                min = feature[i]->size;
          //  if (i == number_of_features )
         //       ++telhoev;
        }

        if (n > n_nodes)
        {
            number_of_features = 0;
        }

/*
// check time
end2 = clock();
cout << "f " << number_of_features << " s " << n << " " << (int)difftime(end1,begin) << "\t"
     << (int)difftime(end2,end1) << "\t" << (int)difftime(end2,begin) << endl;
     cout << "\n";
if (number_of_features < 70)
{
    for (int i = 0;i < number_of_features; ++i)
        cout << feature[i]->size << " ";
}
getchar();*/
    }




    //cout << "time for update " << (int)tot << endl;

	// to include the pulse number 0
	++pulse_index;
	//printf(" hoev %d ",telhoev);
    return pulse_index;
}

int *ReconstructGraph(int n_nodes, PGraphNode *&DPT_Graph, int *pulseRange, int size_pulseRange = 0, int offset = 0)
{
    PGraphNode *CPulse = 0;
    int height = 0;

    int* output = 0;
    output = new int[n_nodes];

    for (int i = 0; i < n_nodes; ++i)
    {
        CPulse = &DPT_Graph[i];

        height = offset;

        while (CPulse->Connected_Pulse != 0)
        {
            CPulse = CPulse->Connected_Pulse;

            for (int boundary = 0; boundary < size_pulseRange; boundary += 2)
            {
                if ((CPulse->size >= pulseRange[boundary])&&(CPulse->size <= pulseRange[boundary + 1]))
                {
                    height += CPulse->height;
                }
            }
            output[i] = height;
        }
    }

    return output;
}

int **PGraph2IntStruct(int n_pulses, PGraphNode *&DPT_Graph,const char* text_output = 0)
{
    PGraphNode *CPulse;
    int size;
    int height;
    int Connected_Pulse_index;
    int Number_of_construction_pulses;

    int **output;
    output = new int*[n_pulses];

    int vector = 3;

    ofstream text_out;

    if (text_output != 0)
        text_out.open(text_output);

    for (int i=0; i < n_pulses; ++i)
    {
        vector = 5;

        CPulse = &(DPT_Graph[(int)i]);
        size = CPulse->size;
        height = CPulse->height;
        if (CPulse->Connected_Pulse != 0)
            Connected_Pulse_index = CPulse->Connected_Pulse->index + 1;
        else
            Connected_Pulse_index = -1;

        Number_of_construction_pulses = CPulse->n_construct;

        // Field 4 - Index of construction pulses
        if (size == 1)
        {
            output[i] = new int[vector];
            output[i][4] = Number_of_construction_pulses + 1; // index of node

            Number_of_construction_pulses = 0;
        }
        else
        {
            vector += Number_of_construction_pulses;
            output[i] = new int[vector];
            output[i][4] = Number_of_construction_pulses;
        }
        // Field 0 - show vector length
        output[i][0] = vector;
        // Field 1 - size
        output[i][1] = size;
        // Field 2 - height
        output[i][2] = height;
        // Field 3 - Connected Pulse
        output[i][3] = Connected_Pulse_index;
        //Field 4 - n -> number of construction pulses or index of pulse if base node
        //Field 5..n - construction pulse index

        for (int indexer = 0; indexer < Number_of_construction_pulses; ++indexer )
        {
            output[i][5 + indexer] = CPulse->ptr_construct[(int)indexer]->index + 1;
        }

        if (text_output != 0)
        {
            for (int g = 0; g < vector - 1; ++g)
                text_out << output[i][g] << "\t";
            text_out << output[i][vector - 1] << "\n";
        }
    }
    if (text_output != 0)
        text_out.close();

    return output;
}

int DPT2PGraph(int *&connectivity, int *&data_structure, PGraphNode *&DPT_Graph)
{
    // timing output
    clock_t begin;
    clock_t end;

    // index of nodes
    WGraphNode **index;
    // index of features
    WGraphNode **feature;
    int number_of_features = 0;
    // datastructure of signal nodes
    WGraphNode *node;
    // the total number of nodes
    int n_nodes = 1;
    int n_pulses = 0;

    // determine number of samples
    for (int k = 0; k < connectivity[0]; ++k)
        n_nodes *= connectivity[k + 1];

    // Create Graph
        // new n_nodes include zero node as well.
        n_nodes = CreateGraph(n_nodes,connectivity,data_structure,index,node,DPT_Graph);
        //n_index = n_nodes;
        number_of_features = OptimizeGraph(n_nodes,feature,index,node);

    // The DPT decomposition
         n_pulses = DPT_decomposition(number_of_features, n_nodes, feature, DPT_Graph,index);

    return n_pulses;
}

int DPT2IntStruct(int *&connectivity, int *&data_structure, int **&DPT_Graph_int, const char *text_output = 0)
{
    PGraphNode *DPT_Graph;

    int n_pulses = DPT2PGraph(connectivity, data_structure, DPT_Graph);

    DPT_Graph_int = PGraph2IntStruct(n_pulses, DPT_Graph,text_output);

    return n_pulses;
}









