/**
 * Copyright (c) 2017 Melown Technologies SE
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * *  Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * *  Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
/** Quadric decimater that computes error in other (physical) system.
 *  Based on OpenMesh::Decimater::ModQuadricT
 */

/*===========================================================================* \
 *                                                                           *
 *                               OpenMesh                                    *
 *      Copyright (C) 2001-2012 by Computer Graphics Group, RWTH Aachen      *
 *                           www.openmesh.org                                *
 *                                                                           *
 *---------------------------------------------------------------------------*
 *  This file is part of OpenMesh.                                           *
 *                                                                           *
 *  OpenMesh is free software: you can redistribute it and/or modify         *
 *  it under the terms of the GNU Lesser General Public License as           *
 *  published by the Free Software Foundation, either version 3 of           *
 *  the License, or (at your option) any later version with the              *
 *  following exceptions:                                                    *
 *                                                                           *
 *  If other files instantiate templates or use macros                       *
 *  or inline functions from this file, or you compile this file and         *
 *  link it with other files to produce an executable, this file does        *
 *  not by itself cause the resulting executable to be covered by the        *
 *  GNU Lesser General Public License. This exception does not however       *
 *  invalidate any other reasons why the executable file might be            *
 *  covered by the GNU Lesser General Public License.                        *
 *                                                                           *
 *  OpenMesh is distributed in the hope that it will be useful,              *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of           *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            *
 *  GNU Lesser General Public License for more details.                      *
 *                                                                           *
 *  You should have received a copy of the GNU LesserGeneral Public          *
 *  License along with OpenMesh.  If not,                                    *
 *  see <http://www.gnu.org/licenses/>.                                      *
 *                                                                           *
\*===========================================================================*/

/*===========================================================================*\
 *                                                                           *
 *   $Revision: 737 $                                                         *
 *   $Date: 2012-10-08 09:33:20 +0200 (Mo, 08 Okt 2012) $                   *
 *                                                                           *
\*===========================================================================*/


//=============================================================================
//
//  CLASS ModQuadricT
//
//=============================================================================

#ifndef geometry_detail_hybrid_decimater_hpp_included_
#define geometry_detail_hybrid_decimater_hpp_included_


//== INCLUDES =================================================================

#include <float.h>
#include <OpenMesh/Tools/Decimater/ModBaseT.hh>
#include <OpenMesh/Core/Utils/Property.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
#include <OpenMesh/Core/Geometry/QuadricT.hh>

#include <boost/optional.hpp>
#include <boost/utility/in_place_factory.hpp>

#include "math/geometry_core.hpp"

#include "dbglog/dbglog.hpp"


//== NAMESPACE ================================================================

namespace geometry { namespace detail {


//== CLASS DEFINITION =========================================================


/** \brief Mesh decimation module computing collapse priority based on error quadrics.
 *
 *  This module can be used as a binary and non-binary module.
 */
template <class MeshT>
class ModQuadricHybrid : public OpenMesh::Decimater::ModBaseT<MeshT>
{
public:

  // Defines the types Self, Handle, Base, Mesh, and CollapseInfo
  // and the memberfunction name()
  DECIMATING_MODULE( ModQuadricHybrid, MeshT, QuadricHybrid );

public:

  /** Constructor
   *  \internal
   */
  ModQuadricHybrid( MeshT &_mesh )
    : Base(_mesh, false)
  {
    unset_max_err();
    Base::mesh().add_property( quadrics_ );
  }


  /// Destructor
  virtual ~ModQuadricHybrid()
  {
    Base::mesh().remove_property(quadrics_);
  }


public: // inherited

  /// Initalize the module and prepare the mesh for decimation.
  virtual void initialize(void);

  /** Compute collapse priority based on error quadrics.
   *
   *  \see ModBaseT::collapse_priority() for return values
   *  \see set_max_err()
   */
  virtual float collapse_priority(const CollapseInfo& _ci)
  {
    using namespace OpenMesh;

    typedef Geometry::QuadricT<double> Q;

    Q q = Base::mesh().property(quadrics_, _ci.v0);
    q += Base::mesh().property(quadrics_, _ci.v1);

    double err = q(vertex(_ci.v1));

    return float( (err < max_err_) ? err : float( Base::ILLEGAL_COLLAPSE ) );
  }


  /// Post-process halfedge collapse (accumulate quadrics)
  virtual void postprocess_collapse(const CollapseInfo& _ci)
  {
    Base::mesh().property(quadrics_, _ci.v1) +=
      Base::mesh().property(quadrics_, _ci.v0);
  }

  /// set the percentage of maximum quadric error
  void set_error_tolerance_factor(double _factor);

  /** Set alternative set of vertices.
   */
  void setAlternativeVertices(const math::Points3d *altVertices) {
      if (!altVertices) {
          altVertices_ = boost::none;
          return;
      }

      altVertices_ = boost::in_place();

      auto &out(*altVertices_);
      for (const auto &v : *altVertices) {
          out.emplace_back(v(0), v(1), v(2));
      }
  }

public: // specific methods

  /** Set maximum quadric error constraint and enable binary mode.
   *  \param _err    Maximum error allowed
   *  \param _binary Let the module work in non-binary mode in spite of the
   *                 enabled constraint.
   *  \see unset_max_err()
   */
  void set_max_err(double _err, bool _binary=true)
  {
    max_err_ = _err;
    Base::set_binary(_binary);
  }

  /// Unset maximum quadric error constraint and restore non-binary mode.
  /// \see set_max_err()
  void unset_max_err(void)
  {
    max_err_ = DBL_MAX;
    Base::set_binary(false);
  }

  /// Return value of max. allowed error.
  double max_err() const { return max_err_; }


private:
  OpenMesh::Vec3d vertex(const typename Mesh::VertexHandle &vh) {
      if (altVertices_) {
          return (*altVertices_)[vh.idx()];
      }
      return OpenMesh::vector_cast<OpenMesh::Vec3d>(Base::mesh().point(vh));
  }

  // maximum quadric error
  double max_err_;

  // this vertex property stores a quadric for each vertex
  OpenMesh::VPropHandleT< OpenMesh::Geometry::QuadricT<double> >  quadrics_;

  typedef std::vector<OpenMesh::Vec3d> VertexList;

  boost::optional<VertexList> altVertices_;
};

//== IMPLEMENTATION ==========================================================


template<class DecimaterType>
void
ModQuadricHybrid<DecimaterType>::
initialize()
{
  using OpenMesh::Geometry::Quadricd;
  // alloc quadrics
  if (!quadrics_.is_valid())
    Base::mesh().add_property( quadrics_ );

  // clear quadrics
  typename Mesh::VertexIter  v_it  = Base::mesh().vertices_begin(),
                             v_end = Base::mesh().vertices_end();

  for (; v_it != v_end; ++v_it)
    Base::mesh().property(quadrics_, v_it).clear();

  // calc (normal weighted) quadric
  typename Mesh::FaceIter          f_it  = Base::mesh().faces_begin(),
                                   f_end = Base::mesh().faces_end();

  typename Mesh::FaceVertexIter    fv_it;
  typename Mesh::VertexHandle      vh0, vh1, vh2;
  typedef OpenMesh::Vec3d          Vec3;
  double                           a,b,c,d, area;

  for (; f_it != f_end; ++f_it)
  {
    fv_it = Base::mesh().fv_iter(f_it.handle());
    vh0 = fv_it.handle();  ++fv_it;
    vh1 = fv_it.handle();  ++fv_it;
    vh2 = fv_it.handle();

    const auto v0(vertex(vh0));
    const auto v1(vertex(vh1));
    const auto v2(vertex(vh2));

    Vec3 n = (v1-v0) % (v2-v0);
    area = n.norm();
    if (area > FLT_MIN)
    {
      n /= area;
      area *= 0.5;
    }

    a = n[0];
    b = n[1];
    c = n[2];
    d = -(v0 | n);

    Quadricd q(a, b, c, d);
    q *= area;

    Base::mesh().property(quadrics_, vh0) += q;
    Base::mesh().property(quadrics_, vh1) += q;
    Base::mesh().property(quadrics_, vh2) += q;
  }
}

//-----------------------------------------------------------------------------

template<class MeshT>
void ModQuadricHybrid<MeshT>::set_error_tolerance_factor(double _factor) {
  if (this->is_binary()) {
    if (_factor >= 0.0 && _factor <= 1.0) {
      // the smaller the factor, the smaller max_err_ gets
      // thus creating a stricter constraint
      // division by error_tolerance_factor_ is for normalization
      double max_err = max_err_ * _factor / this->error_tolerance_factor_;
      set_max_err(max_err);
      this->error_tolerance_factor_ = _factor;

      initialize();
    }
  }
}


template <class MeshT>
class ModQuadricConvexT : public OpenMesh::Decimater::ModQuadricT<MeshT>
{
private:

    float concaveVertexModifier_;

    //  0 - cannot decide
    //  1 - convex
    // -1 - concave
    int check_vertex_convexity (const OpenMesh::VertexHandle &vertex_handle)
    {
        MeshT &omMesh = Base::mesh();
        OpenMesh::PolyConnectivity::ConstVertexVertexCWIter neighbour_it = omMesh.vv_cwiter(vertex_handle);
        if (!neighbour_it.is_valid()) {
            LOG (info3) << "Isolated vertex " << vertex_handle << " skipped.";
            return 0;
        }
        auto start_vertex = omMesh.point(*neighbour_it); // The main point to be used for volume computations
        auto lateral_edge = start_vertex - omMesh.point(vertex_handle);

        auto prev_neighbour_it = ++neighbour_it;
        if (!prev_neighbour_it.is_valid()) {
            LOG (info3) << "Vertex with only one edge " << vertex_handle << " skipped.";
            return 0;
        }
        ++neighbour_it; // Move one vertex forward
        if (!neighbour_it.is_valid()) {
            LOG (info3) << "Vertex with only one face " << vertex_handle << " skipped.";
            return 0;
        }

        // We have the following relation:
        //      start_vertex -> prev_neighbour -> neighbour
        double volume_parallelepiped = 0;
        for ( /*empty*/ ; neighbour_it.is_valid(); ++neighbour_it) {
            auto v0 = omMesh.point (*prev_neighbour_it);
            auto v1 = omMesh.point (*neighbour_it);
            // v0, v1, start_vertex in most of cases do not form a face!
            auto edge1 = start_vertex - v0;
            auto edge2 = start_vertex - v1;

//            LOG (debug) << "(" << lateral_edge << ", [ "<< edge1 << " , " << edge2 <<"])";
            auto a = cross(edge1, edge2);
//            LOG (debug) << "(" << lateral_edge << ", " << a << ")";
            auto b = dot (lateral_edge, a);
//            LOG (debug) << b;

            volume_parallelepiped += b;
        }

        if (volume_parallelepiped > 0 ) {
            //LOG (debug) << "plus -> convex " << volume_parallelepiped/6.0;
            return 1;
        } else {
            //LOG (debug) << "minus -> concave " << volume_parallelepiped/6.0;
            return -1;
        }
    }


public:

  // Defines the types Self, Handle, Base, Mesh, and CollapseInfo
  // and the memberfunction name()
  DECIMATING_MODULE( ModQuadricConvexT, MeshT, QuadricConvex );

public:

  ModQuadricConvexT( MeshT &_mesh)
    : OpenMesh::Decimater::ModQuadricT<MeshT>(_mesh)
  {
      concaveVertexModifier_ = 1; // By default do not do anything
  }

  void set_concave_vertex_modifier(float value) {
      concaveVertexModifier_ = value;
  }

  /// Destructor
  virtual ~ModQuadricConvexT()
  {
  }


public: // inherited

  /**
   *  Compute collapse priority based on error quadrics
   *  and modifies it based on convexity of the vertex
   */
  virtual float collapse_priority(const CollapseInfo& _ci)
  {
    float priority = OpenMesh::Decimater::ModQuadricT<MeshT>::collapse_priority(_ci);
    if (concaveVertexModifier_ == 1) {
        // If we are not going to modify it => do nothing
        return priority;
    }

    // We want to modify priority somehow
    if (priority == Base::ILLEGAL_COLLAPSE) {
        return priority;
    }
    // v0 -> vertex to be removed
    // v1 -> remaining vertex
    auto result = check_vertex_convexity (_ci.v0);
    if (result == -1) {
//        LOG (debug) << "Point " << Base::mesh().point (_ci.v0) << " is concave";
        return priority * concaveVertexModifier_;
    }
    return priority;
  }

};

//=============================================================================
} } // namespace geometry::detail
//=============================================================================

//=============================================================================
#endif // geometry_detail_hybrid_decimater_hpp_included_
//=============================================================================

