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
#include "volumeop.hpp"
#include "math/geometry_core.hpp"

namespace geometry{

//morphological operations
StructuralElement::StructuralElement(SEShape shape, float size){
    switch(shape){
    case SEShape::CUBE:
        {
            int rectHalfSize = std::floor(size/2);
            for(int rz=-rectHalfSize; rz<= rectHalfSize; ++rz){
                for(int ry=-rectHalfSize; ry<= rectHalfSize; ++ry){
                    for(int rx=-rectHalfSize; rx<= rectHalfSize; ++rx){
                        math::Point3i offset(rx,ry,rz);
                        offsets_.push_back(offset);
                    }
                }
            }
        }
        break;
    case SEShape::SPHERE:
        {
            int rectHalfSize = std::floor(size/2);
            for(int rz=-rectHalfSize; rz<= rectHalfSize; ++rz){
                for(int ry=-rectHalfSize; ry<= rectHalfSize; ++ry){
                    for(int rx=-rectHalfSize; rx<= rectHalfSize; ++rx){
                        math::Point3i offset(rx,ry,rz);
                        float distance=ublas::norm_2(offset);
                        if(distance<=size/2) offsets_.push_back(offset);
                    }
                }
            }
        }
        break;
    };
}


} //namespace geometry

