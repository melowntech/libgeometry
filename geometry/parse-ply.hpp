/**
 * Copyright (c) 2021 Melown Technologies SE
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
 * @file geometry/parse-ply.hpp
 * @author Tomas Novak <tomas.novak@melowntech.com>
 *
 * PLY parser that allows loading ply files using simple interface
 *
 * Specify your own `PlyParserBase` child and load any ply you want!
 */

#ifndef geometry_parseply_hpp_included_
#define geometry_parseply_hpp_included_

#include <sstream>
#include <string>

#include <boost/filesystem/path.hpp>
#include <boost/algorithm/string.hpp>

#include "dbglog/dbglog.hpp"

#include "mesh.hpp"

namespace geometry
{

namespace fs = boost::filesystem;
namespace ba = boost::algorithm;

/**
 * Property of element in PLY file (e.g. x, y, z, list of vertex indicies, ...)
 */
struct PlyProperty
{
    /**
     * Property type
     * 
     * Contains whole "list <numerical-type> <numerical-type>" for lists.
     */
    std::string type;

    /**
     * Property name
     */
    std::string name;

    /**
     * Construct a PlyProperty object
     * 
     * @param[in] type property type
     * @param[in] name property name
     */
    PlyProperty(const std::string& type, const std::string& name)
        : type(type),
          name(name)
    { }
};

/**
 * Element in PLY file (e.g. vertex, face, ...)
 */
struct PlyElement
{
    /**
     * Element name
     */
    std::string name;

    /**
     * Number of elements in file
     */
    std::size_t count;

    /**
     * List of element properties
     */
    std::vector<PlyProperty> props;

    /**
     * Construct a PlyElement object with empty properties
     *
     * @param[in] name element name
     * @param[in] count number of elements in file
     */
    PlyElement(const std::string& name, const std::size_t& count)
        : name(name),
          count(count)
    { }
};

/**
 * Virtual class defining interface for PLY parser
 *
 * Defines a set of functions called by the parser when loading a file. Vertex
 * and face elements have a separate callback as that is what will be mostly
 * used.
 */
class PlyParserBase
{
public:
    virtual ~PlyParserBase() {}

    /**
     * Accepts a list of elements with their properties as defined in header
     *
     * Called exactly once, before calling other methods.
     *
     * @param[in] elements vector containing elements of PLY file
     * @param[in] comments vector containing comments in the header
     */
    virtual void loadHeader(const std::vector<PlyElement>& elements,
                            const std::vector<std::string>& comments)
        = 0;

    /**
     * Reads properties corresponding to a vertex element
     *
     * The function reads all properties of a vertex (as specified in header),
     * but leaves the newline character(s) in place.
     *
     * @param[in] f stream to read the properties from
     */
    virtual void addVertex(std::istream& f) = 0;

    /**
     * Reads properties corresponding to a face element
     * 
     * The function reads all properties of a face (as specified in header),
     * but leaves the newline character(s) in place.
     *
     * @param[in] f stream to read the properties from
     */
    virtual void addFace(std::istream& f) = 0;

    /**
     * Reads properties of elements other than "vertex" and "face"
     *
     * Other elements are not common in PLY files, but possible.
     * The function reads all properties of an element (as specified in header),
     * but leaves the newline character(s) in place.
     *
     * @param[in] f stream to read the properties from
     * @param[in] idx line index to infer element type (line after "end_header"
     *                has idx=0)
     */
    void addOtherElement(std::istream& /* f */,
                         const std::size_t /* idx */)
    {
        LOGTHROW(err4, std::runtime_error)
            << "Parsing other elements is not implemented";
    }
};

/**
 * Parses a ply file, calling functions from the parser object
 * 
 * @param[in] parser parser defining callbacks
 * @param[in] path to ply file
 */
void parsePly(PlyParserBase& parser, const fs::path& path);

}  // namespace geometry

#endif // geometry_parseply_hpp_included_
