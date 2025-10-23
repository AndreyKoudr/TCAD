/**
BSD 2-Clause License

Copyright (c) 2025, Andrey Kudryavtsev (andrewkoudr@hotmail.com)
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*******************************************************************************

  Templated CAD

  tsystems.h

  Solutions to systems of linear equations

*******************************************************************************/

#pragma once

#include "tbasics.h"

namespace tcad {

/**
  Simple solutions to linear, over/underdetermined systems with plain C, no SIMD.
*/

/** Transpose matrix, no bound controls is performed */
template <class T> void transpose(T *matrix, int n1, int n2, T *transpose)
{
  for (int i = 0; i < n1; i++)
  {
    for (int j = 0; j < n2; j++)
    {
      transpose[j * n1 + i] = matrix[i * n2 + j];
    }
  }
}

/** Fast multiplication of matrices x3[..,m3] = 
  x1[..,m1] * x2[..,m2]
  where m1,m2,m3 - SECOND dimensions of arrays x1,x2,x3 respectively;
  if array if one-dimensional, then the m should be equal to this 
  dimension, multiplied areas have the following dimensions : n1 x n2 x n3 */
template <class T> void FastM(T *x1, T *x2, T *x3, int n1, int n2, int n3,
  int m1, int m2, int m3)
{
  int m1len,m2len,m3len;
  T *P1,*P2,*P3,t;

  m1len = m1; m2len = m2; m3len = m3;

  for (int i = 0; i < n1; i++)
  {
    P3 = x3; P3 += (i * m3len - 1);
    for (int j = 0; j < n3; j++)
    {
      P3++;
      P1 = x1; P1 += (i * m1len - 1);
      P2 = x2; P2 += (j - m2len); 
      t = 0;
      for (int k = 0; k < n2; k++)
      { 
        P1++;
        P2 += m2len;
        t += (*P1) * (*P2);
      };
      *P3 = t;
    }
  }
}
                
/** Solve square system of equations */
template <class T> bool solveSystem(int N, T A[], T B[], T tolerance)
{
  int K,K1,J,I;
  int I0,I1,I2;
  T AKK;

  for (K = 0; K < N; K++)
  {
    K1 = K + 1;
    AKK = A[K*N+K];

    if (std::abs(AKK) < tolerance)
	  {	
	    return false;
    }

    B[K] /= AKK;
    if (K == (N-1)) break;

    for (J = K1; J < N; J++)
    {
      I0 = K*N+J;
      A[I0] /= AKK;
      for (I = K1; I < N; I++) 
      {
        I1 = I*N+J;
        I2 = I*N+K;
        A[I1] -= (A[I2] * A[I0]);
      }
      B[J] -= (A[J*N+K] * B[K]);
    }
  }

Back:
  K1 = K;
  K -= 1;
  if (K < 0) goto Success;
  for (J = K1; J < N; J++) B[K] -= A[K*N+J] * B[J];
  goto Back;

Success :
  return true;
}

/** Reshuffle matrices A and B from row */
template <class T, class Tint> void findPivot(Tint N, Tint row, T A[], T B[], std::vector<T> &temp)
{
  for (Tint K = row; K < N - 1; K++)
  {
    Tint K1 = K + 1;
                          
    // Find equation with the largest diagonal element
    Tint index = K; T dmax = std::abs(A[K*N+K]);
    for (Tint i = K1; i < N; i++)
    {
      Tint dindex = i * N + K;
      if (std::abs(A[dindex]) > dmax)
      {
        index = i;
        dmax = std::abs(A[dindex]);
      }
    }
                          
    // Make permutation
    if (index != K)
    {
      // Swap matrix raws
      memmove(&temp[0],&A[K * N],sizeof(T) * N);
      memmove(&A[K * N],&A[index * N],sizeof(T) * N);
      memmove(&A[index * N],&temp[0],sizeof(T) * N);

      // Swap right-hand side values
      T btemp = B[K];
      B[K] = B[index];
      B[index] = btemp;
    }
  }
}

/** Solve linear system of equations Ax = B with pivoting */
template <class T, class Tint> bool solveSystemWithPivoting(Tint N, T A[], T B[], T tolerance)
{
  Tint K,K1,J,I;
  Tint I0,I1,I2;
  T AKK;

  // Temporary storage to swap matrix rows
  std::vector<T> temp(N);

  // temp
  std::vector<T> akk;
  T akkmin;

  // Solution itself
  for (K = 0; K < N; K++)
  {
    // Reshuffle matrices A and B
    findPivot<T,Tint>(N,K,A,B,temp);

    K1 = K + 1;
    AKK = A[K*N+K];

    // Temp
    akk.push_back(AKK);
    if (K == 0)
    {
      akkmin = std::abs(AKK);
    } else
    {
      if (std::abs(AKK) < akkmin) 
        akkmin = std::abs(AKK);
    }

    // Pivot is too small
    if (std::abs(AKK) < tolerance)
    {	
      return false;
    }

    B[K] /= AKK;
    if (K == (N-1)) break;

    for (J = K1; J < N; J++)
    {
      I0 = K*N+J;
      A[I0] /= AKK;
      for (I = K1; I < N; I++) 
      {
        I1 = I*N+J;
        I2 = I*N+K;
        A[I1] -= (A[I2] * A[I0]);
      }
      B[J] -= (B[K] * A[J*N+K]);
    }
  }

  // Back substitution
  for (;;)
  {
    K1 = K;
    K -= 1;
    if (K < 0) 
      break;
    for (J = K1; J < N; J++) 
      B[K] -= B[J] * A[K*N+J];
  }

  return true;
}

/** Solve system 2 x 2. */
template <typename T> bool solveSystem2x2(T a11, T a12, T a21, T a22, T b1, T b2, T *x1, T *x2)
{
  T det = a11 * a22 - a12 * a21;
  if (std::abs(det) <= TOLERANCE(T)) 
    return false;
  *x1 = (b1 * a22 - b2 * a12) / det;
  *x2 = (a11 * b2 - a21 * b1) / det;
  return true;
}

/** Solve system of equations for vectors as unknowns. */
template <class T, class TV> bool solveSystemVec(int N, T A[], TV B[], T tolerance)
{
  int K, K1, J, I;
  int I0, I1, I2;
  T AKK;

  for (K = 0; K < N; K++)
  {
    K1 = K + 1;
    AKK = A[K*N + K];

    if (fabs(AKK) < tolerance)
    {
      return false;
    }

    B[K] /= AKK;
    if (K == (N - 1)) 
      break;

    for (J = K1; J < N; J++)
    {
      I0 = K*N + J;
      A[I0] /= AKK;
      for (I = K1; I < N; I++)
      {
        I1 = I*N + J;
        I2 = I*N + K;
        A[I1] -= (A[I2] * A[I0]);
      }
      B[J] -= (B[K] * A[J*N + K]);
    }
  }

Back:
  K1 = K;
  K -= 1;
  if (K < 0) goto Success;
  for (J = K1; J < N; J++) B[K] -= B[J] * A[K*N + J];
  goto Back;

Success:
  return true;
}

/** Solve system of equations with once pivoting only original matrix */
template <class T> bool solveSystemWithPivotingOnce(int N, T A[], T B[], T tolerance)
{
  int K,K1,J,I;
  int I0,I1,I2;
  T AKK;
  T *temp = (T *) malloc(N * sizeof(T));

  if (temp == nullptr) return false;

                              // reshuffle matrices A and B
  for (int K = 0; K < N - 1; K++)
  {
    K1 = K + 1;
                              // find equation with the biggest diagonal element
    int index = K; T dmax = FABS(A[K*N+K]);
    for (int i = K1; i < N; i++)
    {
      int dindex = i * N + K;
      if (FABS(A[dindex]) > dmax)
      {
        index = i;
        dmax = FABS(A[dindex]);
      };
    };
                              // make permutation
    if (index != K)
    {
                              // swap matrix raws
      memmove(temp,&A[K * N],sizeof(T) * N);
      memmove(&A[K * N],&A[index * N],sizeof(T) * N);
      memmove(&A[index * N],temp,sizeof(T) * N);
                              // swap right-hand side values
      T btemp = B[K];
      B[K] = B[index];
      B[index] = btemp;
    };
  };

  for (K = 0; K < N; K++)
  {
    K1 = K + 1;
    AKK = A[K*N+K];

    if (FABS(AKK) < tolerance)
	  {	
      FREE(temp);
	    return false;
    };

    B[K] /= AKK;
    if (K == (N-1)) break;

    for (J = K1; J < N; J++)
    {
      I0 = K*N+J;
      A[I0] /= AKK;
      for (I = K1; I < N; I++) 
      {
        I1 = I*N+J;
        I2 = I*N+K;
        A[I1] -= (A[I2] * A[I0]);
      };
      B[J] -= (A[J*N+K] * B[K]);
    };
  };

Back:
  K1 = K;
  K -= 1;
  if (K < 0) goto Success;
  for (J = K1; J < N; J++) B[K] -= A[K*N+J] * B[J];
  goto Back;

Success :
  FREE(temp);
  return true;
};

/** Solve overdetermined system, N1(number of rows) x N2(columns);
  solution vector is in first N2 values of B */
template <class T> bool SolveSystemOverdetemined(int N1, int N2, T A[], T B[], T tolerance)
{
                              // make transpose
  T *AT = (T *) malloc(N1 * N2 * sizeof(T));
  if (AT == nullptr)
    return false;
  transpose(A,N1,N2,AT);
                              // AT * A
  T *ATA = (T *) malloc(N2 * N2 * sizeof(T));
  if (ATA == nullptr)
  {
    free(AT);
    return false;
  };
  FastM(AT,A,ATA,N2,N1,N2,N1,N2,N2);
                              // AT * B
  T *ATB = (T *) malloc(N2 * sizeof(T));
  if (ATB == nullptr)
  {
    free(AT);
    free(ATA);
    return false;
  };
  FastM(AT,B,ATB,N2,N1,1,N1,1,1);
                              // free AT
  free(AT);
                              // solve system
  bool res = solveSystemWithPivotingOnce(N2,ATA,ATB,tolerance);
                              // copy to first N2 vectors
  memmove(B,ATB,N2 * sizeof(T));
                              // free
  free(ATA);
  free(ATB);

  return res;
}

/**
  Banded matrix (may be non-symmetric). Essentially this is mapping from [i,j],
  i - row, counted from top to bottom, j - column.
  Matrix is stored by row parts within the band. It stores both asymmetric and symmetric parts.
  For example, bandwidth = 5 :
  XXDXX 

  class TB for the right-hand side
*/
template <class T, class TB> class BandedMatrixSimple {
public:
                              // matrix order
  int height = 0;
                              // original half band width, including diagonal
  int bandwidth = 0;

  /** Constructor. */
  BandedMatrixSimple()
  {
  }

  /** Constructor : n - matrix order, bandwidth - matrix bandwidth (including diagonal). */
  BandedMatrixSimple(int n, int pbandwidth)
  {
    height = n;
    bandwidth = pbandwidth;

    buffer.resize(n * bandwidth);
  }
  
  /** Destructor */
  ~BandedMatrixSimple()
  {
  }

	/** Get element index in buffer; i,j are positions in GLOBAL square matrix, each in the range 
    0..(height - 1) */
	inline int elementIndex(int i, int j)
  {
    assert(i >= 0 && i < height);
    assert(j >= 0 && j < height);

    int halfbandwidth = (bandwidth - 1) / 2;
    int index = i * bandwidth + j + halfbandwidth - i;

    return index;
  }

  /** Set matrix element. */
  void setElement(int i, int j, T value)
  {
    buffer[elementIndex(i,j)] = value;
  }

  /** Get matrix element. */
  T getElement(int i, int j)
  {
    return buffer[elementIndex(i,j)];
  }

  /** Get matrix element. */
  T &element(int i, int j)
  {
    return buffer[elementIndex(i,j)];
  }

  /** Element number OK?. */
  bool elementOK(int i, int j)
  {
    int index = elementIndex(i,j);
    return (index >= 0 && index < int(buffer.size()));
  }

  /** Solve linear system. */
  bool solveSystem(TB *B, T zero)
  {
    int K,K1,J,I,Ie,Je;
    T AKK;

    int N = height;
    int MB = (bandwidth - 1) / 2 - 1;

    for (K = 0; K < N; K++)
    {
      K1 = K + 1;
      AKK = buffer[elementIndex(K,K)];
      if (fabs(AKK) < zero)
		  {	
			  return false;
      }

      B[K] /= AKK;
      if (K == (N-1)) 
        break;

		  Je = K1 + MB; if (Je > (N - 1)) Je = N - 1;
      for (J = K1; J <= Je; J++)
      {
        buffer[elementIndex(K,J)] /= AKK; // *this[K,J] /= AKK

			  Ie = K1 + MB; if (Ie > (N - 1)) Ie = N - 1;
        for (I = K1; I <= Ie; I++) 
        {
          buffer[elementIndex(I,J)] -= (buffer[elementIndex(I,K)] * (buffer[elementIndex(K,J)])); // *this[I,J] -= *this[I,K] * *this[K,J]
        };
        B[J] -= (B[K] * buffer[elementIndex(J,K)]);
      }
    }

    for (;;)
    {
      K1 = K;
      K -= 1;
      if (K < 0) break;
      int K2 = K1 + MB;
      if (K2 > (N - 1)) K2 = N - 1;
      for (J = K1; J <= K2; J++) 
      {
        B[K] -= (B[J] * buffer[elementIndex(K,J)]);
      }
    }

    return true;
  }

  /** Change equation to all zeroes and 1.0 at diagonal. */
  void degenerateEquation(int index)
  {
    unsigned char *e = (unsigned char *) &buffer[bandwidth * index];
    memset(e,0,bandwidth * sizeof(T));
    buffer[elementIndex(index,index)] = 1.0;
  }

  /** Max diagonal element. */
  T maxDiagonal()
  {
    T max = -std::numeric_limits<T>::max();

    for (int i = 0; i < height; i++)
    {
      T d = buffer[elementIndex(i,i)];
      max = std::max<T>(max,d);
    }

    return max;
  }

private:
                              // memory to hold matrix
  std::vector<T> buffer;
};


}