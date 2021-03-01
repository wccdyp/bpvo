/*
   This file is part of bpvo.

   bpvo is free software: you can redistribute it and/or modify
   it under the terms of the Lesser GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   bpvo is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   Lesser GNU General Public License for more details.

   You should have received a copy of the Lesser GNU General Public License
   along with bpvo.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * Contributor: halismai@cs.cmu.edu
 */

#include "bpvo/warps.h"
#include <cmath>

namespace bpvo {

Matrix44 HartlyNormalization(const typename detail::warp_traits<RigidBodyWarp>::PointVector& pts)
{
  Point c(Point::Zero());
  for(const auto& p : pts)
    c.noalias() += p;
  c /= (float) pts.size();

  float m = 0.0f;
  for(const auto& p : pts)
    m += (p - c).norm();
  m /= (float) pts.size();

  float s = std::sqrt(3.0) / std::max(m, 1e-6f);

  Matrix44 ret;
  ret.block<3,3>(0,0) = s * Matrix33::Identity();
  ret.block<3,1>(0,3) = -s*c.head<3>();
  ret.block<1,3>(3,0).setZero();
  ret(3,3) = 1.0f;

  return ret;
}

}; // bpvo

