/******************************************************************************
Discrete Pulse Transform Library.
DPT Library v0.1

Copyright (C) 2013, Gene Stoltz
All rights reserved.
******************************************************************************/


#if !defined _DPT_H_
    #define _DPT_H_

// See section "Functions for external use" (bottom of h file)


// Not tested
//#define DATA_RANGE        //add entropy calculation in decomposition

/*******************************************************************************
    Structures
*******************************************************************************/
// Pulse Graph Node Structure
// the base nodes is the nodes showing the position of the pulses in the datstructure
struct PGraphNode {
    int size;                       // size of node
    int height;                     // difference pulse height
    PGraphNode *Connected_Pulse;   // the pulse that this pulse help constructing
    int n_construct;                // number of pulses in construction
                                    // if "size = 1" then n_construct denotes the datapoint position
    PGraphNode **ptr_construct;     // adresses of constructed pulses.
    #ifdef DATA_RANGE
        int *histogram;
    #endif
    //int status;                     //
    //void *CLand;                    // variable in progress - unused
    int index;
};

// Work Graph Node Structure
struct WGraphNode {
    int index;                      // position within the index
    int empty;                      // show if node is empty
    WGraphNode *ptr_merged;         // if node is empty store merged node position
    int size;                       // size of node
    int value;                      // value of node - not variance height
    int n_connected;                // number of connected nodes
    WGraphNode **ptr_connected;     // array of memory pointers of connected nodes
    int n_construct;                // number of pulses this node created
    PGraphNode **ptr_construct;     // array of memory pointers to created pulses
    #ifdef DATA_RANGE
        int *histogram;
    #endif
};


/*******************************************************************************
    Internal functions for DPT.cpp
*******************************************************************************/
int CreateGraph(int n_nodes, int *&connectivity, int *&data_structure, WGraphNode **&index, WGraphNode *&node, PGraphNode *&base);

#ifdef DATA_RANGE
int UpdateHistogram(WGraphNode *&CNode, WGraphNode *&TNode);
#endif

extern inline int UpdateConnectedComp(WGraphNode *&CNode, WGraphNode *&TNode);
extern inline int UpdateConstructionComp(WGraphNode *&CNode, WGraphNode *&TNode);
extern inline bool CheckNodeConnectedness(WGraphNode  *&CNode, int &i_component);
extern inline int CheckFeature(int &n_nodes, int &feature_count,WGraphNode *&CNode, WGraphNode **&feature);
int OptimizeGraph(int n_nodes, WGraphNode **&feature, WGraphNode **&index, WGraphNode *&node);
extern inline int OptimizeWorkGraph(WGraphNode *&CNode);
int AddPulse(int output_value, int &pulse_index,WGraphNode *&CNode,WGraphNode *&TNode,PGraphNode *&pulse );
int DPT_decomposition(int &number_of_features, int &n_nodes, WGraphNode **&feature, PGraphNode *&pulse, WGraphNode **&index);


/******************************************************************************
    Functions for external use
*******************************************************************************/
/* SUMMARY

 ReconstructGraph     - Reconstruct data, using pulses of certain sizes
 PGraph2IntStruct     - Convert PGraph structure to Int structure and text file
 DPT2PGraph           - Perform DPT and get PGraph Structure
 DPT2IntStruct        - Perform DPT and get Int struct
                                                                                */

/*
    ReconstructGraph
        Reconstruct data structure using pulseRange for guidance
    pulseRange -> vector of int of size size_pulseRange
                  must be multiples of 2
                  Even index - define lower boundary
                  Odd index - define upper boundary
                  Index start at 0
    offset ->   gives an offset to reconstructed data
                example usage: to prevent negative data values.
                                                                            */
int *ReconstructGraph(int n_nodes, PGraphNode *&DPT_Graph, int *pulseRange, int size_pulseRange, int offset);

// Convert PGraph structure to integer structure of 2d array, also option to output to text file
int **PGraph2IntStruct(int n_pulses, PGraphNode *&DPT_Graph,const char* text_output);

/*
    DPT2PGraph
    Input:      connectivity <- pointer conectivity vector
                data_structure <- pointer to data structure

    Output:     int n_pulses <- total number of pulses including the base nodes/pulses.

    Require:    DPT_Graph <- pointer to where DPT_Graph must be saved
                                                                                        */
int DPT2PGraph(int *&connectivity, int *&data_structure, PGraphNode *&DPT_Graph);
/*
    DPT2IntStruct
        Same as DPT2PGraph but output InStruct
                                                                                        */
int DPT2IntStruct(int *&connectivity, int *&data_structure, int **&DPT_Graph_int, const char *text_output);

#endif
