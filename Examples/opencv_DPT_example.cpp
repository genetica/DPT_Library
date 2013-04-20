/******************************************************************************
Example opencv_DPT_example.cpp
Discrete Pulse Transform Library.
DPT Library v0.1

Copyright (C) 2013, Gene Stoltz
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

 1. Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
    notice, this list of conditions and the following disclaimer in
    the documentation and/or other materials provided with the
    distribution.
 3. Neither the name of skimage nor the names of its contributors may be
    used to endorse or promote products derived from this software without
    specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
******************************************************************************/

#include <iostream>
#include <fstream>

#include <stdio.h>

#include <time.h>

#include "DPT.h"

#include <opencv2/core/core_c.h>
#include <opencv2/highgui/highgui_c.h>

// Input file name
#define FILENAME "chelsea.jpg"

int CheckDecomposition(PGraphNode *&pulse, IplImage *&img);
int CreateDataStructureOpenCV(int Type, int *&connectivity, int *&data_structure, char *filename);
int * IplImage2struct(IplImage *&img);
IplImage *struct2IplImage(int *&data_struct, int width, int height, int depth = 8);

int main()
{
// Variables
    // timing variables
    clock_t begin;
    clock_t end;

    // DPT input
    int *connectivity;
    int *data_structure;


    // data file name variable
    char data_filename[32];
    int Type;

    // the total number of nodes
    int n_nodes = 0;
        //int n_index;
    int n_pulses;
    // The DPT decomposition variable
    PGraphNode *DPT_Graph;



    //Specific for CreateDataStructureOpenCV function
        Type = 1;                       // Type 0 - for 1-d data
                                        // Type 1 - for image data
        connectivity = new int [1];
        connectivity[0] = 4;            // 4 for 4 connectivity
                                        // 8 for 8 connectivity
        strcpy(data_filename,FILENAME); // Copy defined filename into variable

// Create Data Structure, OpenCV specific
    begin=clock();

        n_nodes = CreateDataStructureOpenCV(Type,connectivity,data_structure,data_filename);

    end=clock();
    printf(" %20s ", "CreateDataStructure");
    printf(" %d \t ms \n",(int)difftime(end,begin));

// DPT decomposition
    begin=clock();

        n_pulses = DPT2PGraph(connectivity, data_structure, DPT_Graph);

    end=clock();
    printf(" %20s ", "DPT Decomposition");
    printf(" %d \t ms \n",(int)difftime(end,begin));


    IplImage *img = cvLoadImage(FILENAME,0);

    int offset = 0;
    int size_pulseRange = 2;
    int *pulseRange;
    pulseRange = new int[size_pulseRange];
    pulseRange[0] = 5000;
    pulseRange[1] = n_nodes;

    int *reconstructedIntData_structure;
// Reconstruct Image with pulses from 5000 to n_nodes
    begin=clock();

    // Get data_structure reconstruction
    reconstructedIntData_structure = ReconstructGraph(n_nodes,DPT_Graph,pulseRange,size_pulseRange,offset);
    // convert data_structure to an IplImage of OPENCV
    IplImage *reconstructed = struct2IplImage(reconstructedIntData_structure,img->width,img->height);

    end=clock();

    printf(" %20s ", "Image Reconstructed");
    printf(" %d \t ms \n",(int)difftime(end,begin));
    // Display reconstructed IplImage
    cvShowImage("Pulses 5000 to  135300",reconstructed);
    cvWaitKey(0);
    cvReleaseImage(&reconstructed);
// Test function to see if DPT decomposition was successfull require OpenCV
    printf("TEST FUNC Error = %d\n",CheckDecomposition(DPT_Graph,img));

    // Clear Memory

    for (int k = n_nodes; k < n_pulses; k++)
    {
        delete [] DPT_Graph[k].ptr_construct;
    }

    delete [] DPT_Graph;
    delete [] connectivity;
    delete [] data_structure;

    return 0;
}





int CheckDecomposition(PGraphNode *&pulse, IplImage *&img)
{
    // Check if decompistion successful.
    PGraphNode *CPulse = 0;

    int error = 0;
    int height = 0;

    int x,y;
    uchar *ptr_img = 0;

    for (y = 0; y < img->height; ++y)
    {
        ptr_img = (uchar *)(img->imageData + y*img->widthStep);
        for (x = 0; x < img->width; ++x)
        {
            // get starting pulse
            CPulse = &pulse[y*img->width + x];
            // reset pulse height
            height = 0;
            // reconstruct pixel value
            while (CPulse->Connected_Pulse != 0)
            {
                CPulse = CPulse->Connected_Pulse;
                height += CPulse->height;
            }
            // check in value equal, if not add an error
            if (abs(ptr_img[x] - height) > 0)
                ++error;
        }
    }
    // stuur error terug, as 0 dan DPT suksesvol
    return error;
}

int CreateDataStructureOpenCV(int Type, int *&connectivity, int *&data_structure, char *filename)
{
    // for text data use type 0. first line must specify the width of the data
    // further on all data must be in one column
    int n_nodes = 0;
    if (Type == 0)
    {
        int width = 0;
        // data file
        std::ifstream data;
        data.open(filename);

        char line[256];


        if (connectivity[0] == 0)
        {
             data.getline(line,256);
             width = atoi(line);
        }

        // count number of data points
        n_nodes = 0;
        while (!data.eof())
        {
            n_nodes++;
            data.getline(line,256);
        }

        // reset text file
        data.seekg(0);
        data.clear();

        data_structure = new int [n_nodes];
        int get_num = 0;
        int i = 0;
        while (!data.eof())
        {
            data >> get_num;
            #ifdef USE_ABS
                data_structure[i] = abs(get_num);
            #else
                data_structure[i] = get_num;
            #endif
            i++;
        }

        if (connectivity[0] == 2)
        {
            delete [] connectivity;
            connectivity = new int [12];
            connectivity[0] = 1;            // one dimension
            connectivity[1] = n_nodes;      // n data points
            connectivity[3] = 2;            // 2 connections
            connectivity[4] =  -1;           // left connection
            connectivity[5] =  1;           // right connection
        }
    }
    else
    if (Type == 1)
    {
        // data file
        IplImage *img = cvLoadImage(filename,0);

        int width = img->width;
        int height = img->height;

        // count number of data points
        n_nodes = width*height;

        data_structure = IplImage2struct(img);

        // use four connectivity
        if (connectivity[0] == 4)
        {
            delete [] connectivity;
            connectivity = new int [12];
            connectivity[0] = 2;            // dimensions
            connectivity[1] = img->width;   // x - dim[1] - size
            connectivity[2] = img->height;  // y - dim[2] - size
            connectivity[3] = 4;            // graph connections
                //right
            connectivity[4] =  1;           // connection[1] change in dim[1]
            connectivity[5] =  0;           // connection[1] change in dim[2]
                //left
            connectivity[6] = -1;           // connection[1] change in dim[1]
            connectivity[7] =  0;           // connection[1] change in dim[2]
                //bottom
            connectivity[8] =  0;           //etc
            connectivity[9] =  1;
                //top
            connectivity[10] =  0;
            connectivity[11] = -1;
        }
        else
        if (connectivity[0] == 8)
        {
            delete [] connectivity;
            connectivity = new int [12];
            connectivity[0] = 2;            // dimensions
            connectivity[1] = img->width;   // x - dim[1] - size
            connectivity[2] = img->height;  // y - dim[2] - size
            connectivity[3] = 8;            // graph connections

            connectivity[4] =  -1;          // left
            connectivity[5] =  0;
            connectivity[6] = -1;           // top left
            connectivity[7] = -1;
            connectivity[8] =  0;           // top
            connectivity[9] =  -1;
            connectivity[10] =  1;          // top right
            connectivity[11] =  -1;

            connectivity[10] =  1;          // right
            connectivity[11] =  0;
            connectivity[12] =  1;          // bottom right
            connectivity[13] =  1;
            connectivity[14] =  0;          // bottom
            connectivity[15] =  1;
            connectivity[16] =  -1;          // bottom left
            connectivity[17] =  1;
        }
    }
  //  else
 //   if (connectivity = )

    return n_nodes;
}

int *IplImage2struct(IplImage *&img)
{
    int *img_struct = 0;

    img_struct = new int[img->width*img->height];

    for (int y = 0; y < img->height; ++y)
    {
        uchar *ptr = (uchar *)(img->imageData + y*img->widthStep);
        for (int x = 0; x < img->width; ++x)
        {
            img_struct[y*img->width + x] = ptr[x];
        }
    }

    return img_struct;
}

IplImage *struct2IplImage(int *&data_struct, int width, int height, int depth)
{
    IplImage *img = 0;
    img = cvCreateImage(cvSize(width,height),depth,1);

    for (int y = 0; y < img->height; ++y)
    {
        uchar *ptr = (uchar *)(img->imageData + y*img->widthStep);
        for (int x = 0; x < img->width; ++x)
        {
            ptr[x] = data_struct[y*img->width + x];
        }
    }

    return img;
}

