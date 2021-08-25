#ifndef UTILS_PAD_HH
#define UTILS_PAD_HH

#include <iostream>
#include <cmath>

#include "../types.hh"

#define MAX_VALUE 0
#define MIN_VALUE 1
#define AVG_VALUE 2
#define CONST_VAL 3

namespace padding {

template <typename T>
T* fill_2d_block(T* array, const size_t* dims, size_t blkSize, size_t b0, size_t b1)
{
    T* block = (T*)malloc(sizeof(T) * blkSize * blkSize);

    size_t _idx0 = b0 * blkSize;
    size_t _idx1 = b1 * blkSize;

    for (size_t j = 0; j < blkSize; j++)
    {
        for (size_t i = 0; i < blkSize; i++)
        {
            size_t gi0 = _idx0 + i;
            size_t gi1 = _idx1 + j;
            if (gi1 >= dims[DIM1] or gi0 >= dims[DIM0]) continue;
            block[i + j*blkSize] = array[gi0 + gi1*dims[DIM0]];
        }
    }

    return block;
}

template <typename T>
T* fill_3d_block(T* array, const size_t* dims, size_t blkSize, size_t b0, size_t b1, size_t b2)
{
    T* block = (T*)malloc(sizeof(T) * blkSize * blkSize * blkSize);

    size_t _idx0 = b0 * blkSize;
    size_t _idx1 = b1 * blkSize;
    size_t _idx2 = b2 * blkSize;

    for (size_t k = 0; k < blkSize; k++)
    {
        for (size_t j = 0; j < blkSize; j++)
        {
            for (size_t i = 0; i < blkSize; i++)
            {
                size_t gi0 = _idx0 + i;
                size_t gi1 = _idx1 + j;
                size_t gi2 = _idx2 + k;
                if (gi2 >= dims[DIM2] or gi1 >= dims[DIM1] or gi0 >= dims[DIM0]) continue;
                block[i + j*blkSize + k*blkSize*blkSize] = array[gi0 + gi1*dims[DIM0] + gi2*dims[DIM0]*dims[DIM1]];
            }
        }
    }
return block;
}

template <typename T>
T find_max(T* array, size_t* dims)
{
    T max = array[0];
    for (size_t k = 0; k < dims[DIM2]; k++)
    {
        for (size_t j = 0; j < dims[DIM1]; j++)
        {
            for (size_t i = 0; i < dims[DIM0]; i++)
            {
                max = (array[i + j*dims[DIM0] + k*dims[DIM1]*dims[DIM0]] > max) ? array[i + j*dims[DIM0] + k*dims[DIM1]*dims[DIM0]] : max;
            }
        }
    }

    return max;
}

template <typename T>
T find_min(T* array, size_t* dims)
{
    T min = array[0];
    for (size_t k = 0; k < dims[DIM2]; k++)
    {
        for (size_t j = 0; j < dims[DIM1]; j++)
        {
            for (size_t i = 0; i < dims[DIM0]; i++)
            {
                min = (array[i + j*dims[DIM0] + k*dims[DIM1]*dims[DIM0]] < min) ? array[i + j*dims[DIM0] + k*dims[DIM1]*dims[DIM0]] : min;
            }
        }
    }

    return min;
}

template <typename T>
T find_avg(T* array, size_t* dims)
{
    T running_sum = 0;
    for (size_t k = 0; k < dims[DIM2]; k++)
    {
        for (size_t j = 0; j < dims[DIM1]; j++)
        {
            for (size_t i = 0; i < dims[DIM0]; i++)
            {
                running_sum += array[i + j*dims[DIM0] + k*dims[DIM0]*dims[DIM1]];
            }
        }
    }

    T avg = running_sum / (float)(dims[DIM2]*dims[DIM1]*dims[DIM0]);

    return avg;
}

template <typename T>
T find_pad_value(T* array, int scalarType, size_t* dims, T constValue)
{
    T scalar = 0;
    if (scalarType == MAX_VALUE)
    {
        scalar = find_max<T>(array, dims);
    }
    else if (scalarType == MIN_VALUE)
    {
        scalar = find_min<T>(array, dims);
    }
    else if (scalarType == AVG_VALUE)
    {
        scalar = find_avg<T>(array, dims);
    }
    else if (scalarType == CONST_VAL)
    {
        scalar = constValue;
    }

    return scalar;
}

template <typename T>
T* find_edge_pad_value(T* block, int padType, int nDims, size_t blkSize, T constValue)
{
    auto padValues = new T[nDims];

    if (padType == CONST_VAL)
    {
        for (int i = 0; i < nDims; i++) padValues[i] = constValue;
        return padValues;
    }
    else
    {
        //find the edge of an array
        if (nDims == 1)
        {
            padValues[0] = block[0];
        }
        else if (nDims == 2)
        {
            size_t dims[3] {blkSize, 1, 1};

            T top_edge[blkSize], *topEdgePtr = top_edge;
            T left_edge[blkSize], *leftEdgePtr = left_edge;

            for (size_t i = 0; i < blkSize; i++)
            {
                top_edge[i]  = block[i];
                left_edge[i] = block[i * blkSize];
            }

            padValues[0] = find_pad_value<T>(topEdgePtr, padType, dims, constValue);
            padValues[1] = find_pad_value<T>(leftEdgePtr, padType, dims, constValue);
        }
        else if (nDims == 3)
        {
            size_t dims[3] {blkSize, blkSize, 1};

            T top_edge[blkSize * blkSize], *topEdgePtr;
            T left_edge[blkSize * blkSize], *leftEdgePtr;
            T front_edge[blkSize * blkSize], *frontEdgePtr;
            topEdgePtr = top_edge;
            leftEdgePtr = left_edge;
            frontEdgePtr = front_edge;

            for (size_t j = 0; j < blkSize; j++)
            {
                for (size_t i = 0; i < blkSize; i++)
                {
                    front_edge[i + j * blkSize] = block[i + j * blkSize];
                    left_edge[i + j * blkSize]  = block[i * blkSize + j * blkSize * blkSize];
                    top_edge[i + j * blkSize]   = block[i + j * blkSize * blkSize];
                }
            }

            padValues[0] = find_pad_value<T>(frontEdgePtr, padType, dims, constValue);
            padValues[1] = find_pad_value<T>(leftEdgePtr, padType, dims, constValue);
            padValues[2] = find_pad_value<T>(topEdgePtr, padType, dims, constValue);
        }

        return padValues;
    }
}

template <typename T>
T* block_pad(T* block, int nDims, int padType, size_t blkSize, T constValue, double padValueModifier, T* padValue)
{
    size_t blockLength = pow((blkSize + 1), nDims);
    T* paddedBlock = new T[blockLength];

    if (nDims == 2)
    {
        size_t dims[3] {blkSize, blkSize, 1};

        T padValueLocal = find_pad_value<T>(block, padType, dims, constValue);
        padValueLocal *= padValueModifier; // adjust padding value for prequantization
        padValueLocal = round(padValueLocal);

        for (size_t j = 0; j < blkSize + 1; j++)
        {
            for (size_t i = 0; i < blkSize + 1; i++)
            {
                paddedBlock[i + j * (blkSize + 1)] = (i == 0 || j == 0) ? padValueLocal : 0;
            }
        }
        *padValue = padValueLocal;
    }
    else if (nDims == 3)
    {
        size_t dims[3] {blkSize, blkSize, blkSize};

        T padValueLocal = find_pad_value<T>(block, padType, dims, constValue);
        padValueLocal *= padValueModifier; // adjust padding value for prequantization
        padValueLocal = round(padValueLocal);

        for (size_t k = 0; k < blkSize + 1; k++)
        {
            for (size_t j = 0; j < blkSize + 1; j++)
            {
                for (size_t i = 0; i < blkSize + 1; i++)
                {
                    paddedBlock[i + j * (blkSize + 1) + k * (blkSize + 1) * (blkSize + 1)] = (i == 0 || j == 0 || k == 0) ? padValueLocal : 0;
                }
            }
        }
        *padValue = padValueLocal;
    }

    return paddedBlock;
}

template <typename T>
T* edge_pad(T* block, int nDims, int padType, size_t* dims, size_t blkSize, T constValue, double padValueModifier, T** padValues)
{
    T* padValuesLocal = find_edge_pad_value(block, padType, nDims, blkSize, constValue);
    for (int i = 0; i < nDims; i++)
    {
        padValuesLocal[i] *= padValueModifier; // adjust padding value for prequantization
        padValuesLocal[i] = round(padValuesLocal[i]);
    }
    //T* padBlockPtr;
    size_t blockLength = pow((blkSize + 1), nDims);
    auto paddedBlock = new T[blockLength] {0};

    if (nDims == 2)
    {
        //T paddedBlock[(blkSize + 1) * (blkSize + 1)] {0};
        for (size_t j = 0; j < blkSize + 1; j++)
        {
            for (size_t i = 0; i < blkSize + 1; i++)
            {
                paddedBlock[i + j * (blkSize + 1)] = (i == 0) ? padValuesLocal[0] : (j == 0) ? padValuesLocal[1] : 0;
            }
        }

      //  padBlockPtr = paddedBlock;
    }
    else if (nDims == 3)
    {
        //T paddedBlock[(blkSize + 1) * (blkSize + 1) * (blkSize + 1)] {0};
        for (size_t k = 0; k < blkSize + 1; k++)
        {
            for (size_t j = 0; j < blkSize + 1; j++)
            {
                for (size_t i = 0; i < blkSize + 1; i++)
                {
                    paddedBlock[i + j * (blkSize + 1) + k * (blkSize + 1) * (blkSize + 1)] = (k == 0) ? padValuesLocal[0] : (j == 0) ? padValuesLocal[1] : (i == 0) ? padValuesLocal[2] : 0;
                }
            }
        }

     //   padBlockPtr = paddedBlock;
    }

    *padValues = padValuesLocal;
    return paddedBlock;
}

template <typename T>
T* x_block_pad(int nDims, size_t blkSize, T padValueLocal)
{
//    T* padBlockPtr;
    size_t blockLength = pow((blkSize + 1), nDims);
    auto paddedBlock = new T[blockLength] {0};

    if (nDims == 2)
    {
        //T paddedBlock[(blkSize + 1) * (blkSize + 1)];
        for (size_t j = 0; j < blkSize + 1; j++)
        {
            for (size_t i = 0; i < blkSize + 1; i++)
            {
                paddedBlock[i + j * (blkSize + 1)] = (i == 0 || j == 0) ? padValueLocal : 0;
            }
        }
 //       padBlockPtr = paddedBlock;
    }
    else if (nDims == 3)
    {
        //T paddedBlock[(blkSize + 1) * (blkSize + 1) * (blkSize + 1)];
        for (size_t k = 0; k < blkSize + 1; k++)
        {
            for (size_t j = 0; j < blkSize + 1; j++)
            {
                for (size_t i = 0; i < blkSize + 1; i++)
                {
                    paddedBlock[i + j * (blkSize + 1) + k * (blkSize + 1) * (blkSize + 1)] = (i == 0 || j == 0 || k == 0) ? padValueLocal : 0;
                }
            }
        }
  //      padBlockPtr = paddedBlock;
    }

    return paddedBlock;
    //return padBlockPtr;
}

template <typename T>
T* x_edge_pad(int nDims, size_t blkSize, T* padValues)
{
    size_t blockLength = pow((blkSize + 1), nDims);
    auto paddedBlock = new T[blockLength] {0};

    if (nDims == 2)
    {
        for (size_t j = 0; j < blkSize + 1; j++)
        {
            for (size_t i = 0; i < blkSize + 1; i++)
            {
                paddedBlock[i + j * (blkSize + 1)] = (i == 0) ? padValues[0] : (j == 0) ? padValues[1] : 0;
            }
        }
    }
    else if (nDims == 3)
    {
        for (size_t k = 0; k < blkSize + 1; k++)
        {
            for (size_t j = 0; j < blkSize + 1; j++)
            {
                for (size_t i = 0; i < blkSize + 1; i++)
                {
                    paddedBlock[i + j * (blkSize + 1) + k * (blkSize + 1) * (blkSize + 1)] = (k == 0) ? padValues[0] : (j == 0) ? padValues[1] : (i == 0) ? padValues[2] : 0;
                }
            }
        }

    }

    return paddedBlock;
}


} // namespace padding

#endif
