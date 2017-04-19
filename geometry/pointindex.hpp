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
/**
 * @file pointindex.hpp
 * @author Jakub Cerveny <jakub.cerveny@ext.citationtech.net>
 *
 * Helper class for removing duplicate vertices.
 */

#ifndef geometry_pointindex_hpp_included_
#define geometry_pointindex_hpp_included_

#include <map>

namespace geometry {

//! Helper class for removing duplicate vertices.
//!
template<typename Point>
class PointIndex
{
    std::map<Point, unsigned> map_;
    unsigned next_;

public:
    PointIndex() : next_(0) {}

    //! Assigns (and returns) a unique index for the given point.
    //! Idices are reused for points that have been seen before.
    unsigned assign(const Point& pt)
    {
        auto pair = map_.insert(std::make_pair(pt, next_));
        if (pair.second) next_++;
        return pair.first->second;
    }
};

} // namespace geometry

#endif // geometry_pointindex_hpp_included_
