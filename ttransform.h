/*
BSD 2-Clause License

Copyright (c) 2020, Andrey Kudryavtsev (andrewkoudr@hotmail.com)
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

  ttransform.h

  TTransform

*******************************************************************************/

#pragma once

#include "tpoint.h"
#include "tplane.h"
#include "tmatrix.h"

namespace tcad {

/** TTransform : transformation of coordinates. */
template <typename T> class TTransform { 
public:
  typedef Matrix<T,4,4> matrix4x4;

  /** Transformation matrices. */
  matrix4x4 ForTransform;			      
  matrix4x4 BackTransform;
  bool BackTransformReady;

  /** Which transformations have been applied? */
  bool TranslationMade;
  TPoint<T> Position;

  bool RotationMade;
  TPoint<T> Angles;
	T RotAngle;

  bool ResizeMade;
  TPoint<T> Coefs;

  bool PerspectiveMade;
  T D;

  bool ScreenMade;
  T PixWidth,PixHeight,Width,Height;

  /** Constructors */
  TTransform()
  {
    ForTransform = IdentityMatrix();
    BackTransformReady = false;

    // which transformations have been applied
    TranslationMade = RotationMade = ResizeMade = PerspectiveMade = ScreenMade = false;
    Position = Angles = TPoint<T>(0.0,0.0,0.0);
    Coefs = TPoint<T>(1.0,1.0,1.0);
    D = 0.0;
    RotAngle = 0.0;
    PixWidth = PixHeight = Width = Height = 0.0;
  }

  /** Make a copy */
  void Copy(TTransform<T> *copy);

  /** Apply transformation to point. */
  TPoint<T> applyTransform(TPoint<T> Point);

  /** Apply back transformation to point. */
  TPoint<T> applyBackTransform(TPoint<T> Point);

  /** Set identity matrix. */                               
  void LoadIdentity();

  /** THE FUNCTIONS BELOW CAN BE CALLED ONLY ONCE TO DEFINE TRANSFORMATION */
  void Translate(TPoint<T> Position);

  /** Set rotation
    angles.X - different
          .Y - course
          .Z - roll */
  void Rotate(TPoint<T> Angles);

  /** Rotate, angle in radians. */
  void Rotate(TPoint<T> U, T Angle);

  /** Resize. */
  void Resize(TPoint<T> Coefs);

  /** Perspective. */
  void Perspective(T d);

  /** Shadow. */
  void Shadow(TPlane<T> *plane, TPoint<T> *dir);

  /** Pixel screen cooridinates in (Z/W,Y/W); */
  void Screen(T pixwidth, T pixheight, T width, T height);

  /** Prepare back transformation matrix. */
  void PrepareBackTransform();

  /** Reverse transform. Call only once! - it replaces tranform matrix by reverse. */
  void Reverse();

  /** Only rotation is made. */
  bool RotationOnly();

  /** Add matrix from another transform. */
  void Add(TTransform &t);

private:

  /** Calculate transformation matrices */
  matrix4x4 IdentityMatrix()
	{
		matrix4x4 R; R.setIdentity();
		return R;
	}

  matrix4x4 TranslationMatrix(TPoint<T> Position)
	{
		matrix4x4 R; R.setTranslation(Position.X,Position.Y,Position.Z);
		return R;
	}

  matrix4x4 RotationMatrix(TPoint<T> Angles)
	{
		matrix4x4 R; R.setRotation(Angles);
		return R;
	}

  matrix4x4 ShadowMatrix(TPlane<T> *plane, TPoint<T> *dir)
	{
		matrix4x4 R; R.setShadow(plane,dir);
		return R;
	}

  matrix4x4 RotationMatrix(TPoint<T> Axis, T Angle)
	{
		matrix4x4 R; 
		if (!Axis > 0.00000001)
		{
			R.setRotation(Axis.X,Axis.Y,Axis.Z,Angle);
		} else
		{
			R.setIdentity();
		}
		return R;
	}

  matrix4x4 ResizeMatrix(TPoint<T> Coefs)
	{
		matrix4x4 R; R.setScaling(Coefs.X,Coefs.Y,Coefs.Z);
		return R;
	}

  matrix4x4 ScreenMatrix(T pixwidth, T pixheight, 
    T width, T height)
	{
		matrix4x4 R; R.setScalingTranslation(1,- pixheight / height,
      pixwidth / width,0,pixheight * 0.5,pixwidth * 0.5);
		return R;
	}

  matrix4x4 PerspectiveMatrix(T d)
	{
		matrix4x4 R; R.setPerspective(d);
		return R;
	}
};

template <typename T> TPoint<T> TTransform<T>::applyTransform(TPoint<T> Point)
{
  Point.W = 1.0;
  return ForTransform * Point;
}

template <typename T> TPoint<T> TTransform<T>::applyBackTransform(TPoint<T> Point)
{
                               // calculate back transformation
  if (!BackTransformReady) PrepareBackTransform();

  Point.W = 1.0;
  return BackTransform * Point;
}

template <typename T> void TTransform<T>::PrepareBackTransform()
{
  BackTransform = (+ForTransform);
  BackTransformReady = true;
}

template <typename T> void TTransform<T>::Reverse()
{
  PrepareBackTransform();
  ForTransform = BackTransform;
}

template <typename T> void TTransform<T>::LoadIdentity()
{
  ForTransform = IdentityMatrix();

  TranslationMade = RotationMade = ResizeMade = false;
  Position = Angles = TPoint<T>(0.0,0.0,0.0);
  Coefs = TPoint<T>(1.0,1.0,1.0);
  BackTransformReady = false;
}

template <typename T> void TTransform<T>::Translate(TPoint<T> position)
{
  matrix4x4 R = TranslationMatrix(position);
  ForTransform = R * ForTransform;

  TranslationMade = true;
  Position = position;
  BackTransformReady = false;
}

template <typename T> void TTransform<T>::Shadow(TPlane<T> *plane, TPoint<T> *dir)
{
  matrix4x4 R = ShadowMatrix(plane,dir);
  ForTransform = R * ForTransform;
  BackTransformReady = false;
}

template <typename T> void TTransform<T>::Rotate(TPoint<T> angles)
{
  matrix4x4 R = RotationMatrix(angles);
  ForTransform = R * ForTransform;

  RotationMade = true;
  Angles = angles;
  BackTransformReady = false;
}

template <typename T> void TTransform<T>::Rotate(TPoint<T> U, T Angle)
{
  matrix4x4 R = RotationMatrix(U,Angle);
  ForTransform = R * ForTransform;

  RotationMade = true;
  Angles = U;
	RotAngle = Angle;
  BackTransformReady = false;
}

template <typename T> void TTransform<T>::Resize(TPoint<T> coefs)
{
  matrix4x4 R = ResizeMatrix(coefs);
  ForTransform = R * ForTransform;

  ResizeMade = true;
  Coefs = coefs;
  BackTransformReady = false;
}

template <typename T> void TTransform<T>::Perspective(T d)
{
  matrix4x4 R = PerspectiveMatrix(d);
  ForTransform = R * ForTransform;

  PerspectiveMade = true;
  D = d;
  BackTransformReady = false;
}

template <typename T> void TTransform<T>::Screen(T pixwidth, T pixheight, T width, T height)
{
  matrix4x4 R = ScreenMatrix(pixwidth,pixheight,width,height);
  ForTransform = R * ForTransform;

  ScreenMade = true;
  PixWidth = pixwidth; PixHeight = pixheight; Width = width; Height = height;
  BackTransformReady = false;
}

#define COPY(v) v = copy->v 

template <typename T> void TTransform<T>::Copy(TTransform<T> *copy)
{
                              // transformation matrices
  COPY(ForTransform);
  COPY(BackTransform);
  COPY(BackTransformReady);

                              // which transformations have been applied
  COPY(TranslationMade);
  COPY(Position);

  COPY(RotationMade);
  COPY(Angles);
	COPY(RotAngle);

  COPY(ResizeMade);
  COPY(Coefs);

  COPY(PerspectiveMade);
  COPY(D);

  COPY(ScreenMade);
  COPY(PixWidth);
  COPY(PixHeight);
  COPY(Width);
  COPY(Height);
}

template <typename T> bool TTransform<T>::RotationOnly()
{
  return (!TranslationMade && RotationMade && !ResizeMade && !PerspectiveMade && !ScreenMade);
}

template <typename T> void TTransform<T>::Add(TTransform<T> &t)
{
  ForTransform = t.ForTransform * ForTransform;

  BackTransformReady = false;
}

#undef COPY

}