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

