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
 *  @file geometry/detail/parse-obj.hpp
 *  @author Vaclav Blazek <vaclav.blazek@citationtech.net>
 *
 *  Boost.Spirit-based OBJ file format parser (implementation).
 */

#ifndef geometry_detail_objparser_hpp_included_
#define geometry_detail_objparser_hpp_included_

#include <boost/config/warning_disable.hpp>
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_object.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/spirit/include/phoenix_fusion.hpp>
#include <boost/spirit/include/phoenix_container.hpp>
#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/spirit/include/qi_match.hpp>
#include <boost/spirit/include/qi_match_auto.hpp>
#include <boost/spirit/include/qi_stream.hpp>

#include <geometry/parse-obj.hpp>

namespace geometry { namespace detail {

struct Vertex : ObjParserBase::Vector3d {};
struct Normal : ObjParserBase::Vector3d {};
struct Texture : ObjParserBase::Vector3d {};
struct MaterialLibrary : std::string {};
struct UseMaterial : std::string {};

typedef ObjParserBase::Facet Facet;

class Obj {
public:
    Obj(ObjParserBase &p)
        : p_(&p), vCount_(), tCount_(), nCount_()
    {}

    Obj& operator+=(const detail::Vertex &v) {
        ++vCount_;
        p_->addVertex(v);
        return *this;
    }

    Obj& operator+=(const detail::Texture &t) {
        ++tCount_;
        p_->addTexture(t);
        return *this;
    }

    Obj& operator+=(const detail::Normal &n) {
        ++nCount_;
        p_->addNormal(n);
        return *this;
    }

    Obj& operator+=(Facet f) {
        // convert to zero-based indices, also support negative indices
        for (auto &v : f.v) {
            v = (v < 0) ? (vCount_ + v) : (v - 1);
        }

        for (auto &t : f.t) {
            t = (t < 0) ? (tCount_ + t) : (t - 1);
        }

        for (auto &n : f.n) {
            n = (n < 0) ? (nCount_ + n) : (n - 1);
        }

        p_->addFacet(f);
        return *this;
    }

    Obj& operator+=(const detail::MaterialLibrary &l) {
        p_->materialLibrary(l);
        return *this;
    }

    Obj& operator+=(const detail::UseMaterial &m) {
        p_->useMaterial(m);
        return *this;
    }

private:
    ObjParserBase *p_;
    int vCount_;
    int tCount_;
    int nCount_;
};

} } // namespace geometry::detail

BOOST_FUSION_ADAPT_STRUCT(
    geometry::detail::Vertex,
    (double, x)
    (double, y)
    (double, z)
)

BOOST_FUSION_ADAPT_STRUCT(
    geometry::detail::Texture,
    (double, x)
    (double, y)
    (double, z)
)

BOOST_FUSION_ADAPT_STRUCT(
    geometry::detail::Normal,
    (double, x)
    (double, y)
    (double, z)
)

BOOST_FUSION_ADAPT_STRUCT(
    geometry::detail::Facet,
    (int, v[0])(int, t[0])(int, n[0])
    (int, v[1])(int, t[1])(int, n[1])
    (int, v[2])(int, t[2])(int, n[2])
)

namespace geometry { namespace detail {

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;

template <typename Iterator, typename Skipper>
struct vertex_parser : qi::grammar<Iterator, Vertex(), Skipper>
{
    vertex_parser() : vertex_parser::base_type(start)  {
        using qi::auto_;

        start %= 'v' >> auto_ >> auto_ >> auto_;
    }

    qi::rule<Iterator, Vertex(), Skipper> start;
};

template <typename Iterator, typename Skipper>
struct texture_parser : qi::grammar<Iterator, Texture(), Skipper>
{
    texture_parser() : texture_parser::base_type(start)  {
        using qi::lit;
        using qi::auto_;

        start %= lit("vt") >> auto_ >> -auto_ >> -auto_;
    }

    qi::rule<Iterator, Texture(), Skipper> start;
};

template <typename Iterator, typename Skipper>
struct normal_parser : qi::grammar<Iterator, Normal(), Skipper>
{
    normal_parser() : normal_parser::base_type(start)  {
        using qi::lit;
        using qi::auto_;

        start %= lit("vn") >> auto_ >> auto_ >> auto_;
    }

    qi::rule<Iterator, Normal(), Skipper> start;
};

template <typename Iterator, typename Skipper>
struct facet_parser : qi::grammar<Iterator, Facet(), Skipper>
{
    facet_parser() : facet_parser::base_type(start)  {
        using qi::auto_;
        using qi::omit;
        using qi::no_skip;
        using qi::lexeme;

        start %= 'f'
            >> lexeme[auto_ >> -('/' >> -auto_) >> -('/' >> -auto_)
                      >> omit[+ascii::space]
                      >> auto_ >> -('/' >> -auto_) >> -('/' >> -auto_)
                      >> omit[+ascii::space]
                      >> auto_ >> -('/' >> -auto_) >> -('/' >> -auto_)
                      ];
    }

    qi::rule<Iterator, Facet(), Skipper> start;
};

template <typename Iterator, typename Skipper>
struct materialLibrary_parser
    : qi::grammar<Iterator, MaterialLibrary(), Skipper>
{
    materialLibrary_parser() : materialLibrary_parser::base_type(start)  {
        using qi::char_;
        using qi::lexeme;

        start %= "mtllib" >> lexeme[+(char_ - ascii::space)];
    }

    qi::rule<Iterator, MaterialLibrary(), Skipper> start;
};

template <typename Iterator, typename Skipper>
struct useMaterial_parser
    : qi::grammar<Iterator, UseMaterial(), Skipper>
{
    useMaterial_parser() : useMaterial_parser::base_type(start)  {
        using qi::char_;
        using qi::no_skip;

        start %= "usemtl" >> qi::lexeme[+(char_ - ascii::space)];
    }

    qi::rule<Iterator, UseMaterial(), Skipper> start;
};

template <typename Iterator, typename Skipper>
struct Obj_parser : qi::grammar<Iterator, Obj(), Skipper>
{
    Obj_parser() : Obj_parser::base_type(start)  {
        using qi::omit;

        vertex %= vertex_parser<Iterator, Skipper>();
        texture %= texture_parser<Iterator, Skipper>();
        normal %= normal_parser<Iterator, Skipper>();
        facet %= facet_parser<Iterator, Skipper>();
        materialLibrary %= materialLibrary_parser<Iterator, Skipper>();
        useMaterial %= useMaterial_parser<Iterator, Skipper>();

        start %= omit[*(vertex[qi::_val += qi::_1]
                        | texture[qi::_val += qi::_1]
                        | normal[qi::_val += qi::_1]
                        | facet[qi::_val += qi::_1]
                        | materialLibrary[qi::_val += qi::_1]
                        | useMaterial[qi::_val += qi::_1]
                        )]
            [qi::_val];
    }

    qi::rule<Iterator, Obj(), Skipper> start;

    vertex_parser<Iterator, Skipper> vertex;
    texture_parser<Iterator, Skipper> texture;
    normal_parser<Iterator, Skipper> normal;
    facet_parser<Iterator, Skipper> facet;
    materialLibrary_parser<Iterator, Skipper> materialLibrary;
    useMaterial_parser<Iterator, Skipper> useMaterial;
};

template <typename Iterator>
struct skipper : qi::grammar<Iterator>
{
    typedef skipper<Iterator> type;

    skipper() : skipper::base_type(start)  {
        comment %= '#' >> *(qi::char_ - qi::char_("\r\n"))
                       >> (qi::eol | qi::eoi);

        start %= ascii::space | comment;
    }

    qi::rule<Iterator> start;
    qi::rule<Iterator> comment;
};

template <typename Iterator>
bool parse(Iterator begin, Iterator end, ObjParserBase &p)
{
    typedef detail::skipper<Iterator> skipper_type;

    detail::Obj_parser<Iterator, skipper_type> qrammar;
    Obj o(p);

    bool r(phrase_parse(begin, end, qrammar, skipper_type(), o));

    return (r && (begin == end));
}

template<typename CharT, typename Traits>
bool parse(std::basic_istream<CharT, Traits> &is, ObjParserBase &p)
{
    is.unsetf(std::ios::skipws);
    typedef boost::spirit::istream_iterator iterator_type;
    return parse(iterator_type(is), iterator_type(), p);
}

} } // namespace geometry::detail


#endif // geometry_detail_objparser_hpp_included_
