/**
 *  @file geometry/parse-obj.cpp
 *  @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 *  Boost.Spirit-based OBJ file format parser.
 */

#include <cstdio>

#include "utility/expect.hpp"

#include "./detail/parse-obj.hpp"

namespace geometry {

bool ObjParserBase::parse(std::istream &is)
{
    return detail::parse(is, *this);
}

bool ObjParserBase::parse(const boost::filesystem::path &path)
{
    LOG(info1) << "Loading OBJ file from " << path << ".";

    std::ifstream f(path.string());
    if (!f.good()) {
        return false;
    }

    f.exceptions(std::ios::badbit | std::ios::failbit);
    auto res(parse(f));
    f.close();
    return res;
}

} // namespace geometry
