/**
 *  @file geometry/parse-obj.cpp
 *  @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 *  Boost.Spirit-based OBJ file format parser.
 */

#include <geometry/detail/parse-obj.hpp>

namespace geometry {

bool ObjParserBase::parse(std::istream &is)
{
    return detail::parse(is, *this);
}

} // namespace geometry
