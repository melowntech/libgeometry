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
 * @file geometry/parse-ply.cpp
 * @author Tomas Novak <tomas.novak@melowntech.com>
 *
 * PLY parser that allows loading ply files using simple interface
 * 
 * Specify your own `PlyParserBase` child and load any ply you want!
 */

#include "parse-ply.hpp"

namespace geometry
{

/// Header parsing functions

/// Helper that throws and exception if EOF bit is set
inline void throwOnEof(const std::ifstream& f)
{
    if (f.eof())
    {
        LOGTHROW(err4, std::runtime_error)
            << "Unexpected EOF reached while parsing PLY header";
    }
}

/// Splits property line to name and type
PlyProperty parseProperty(const std::string& line)
{
    PlyProperty prop;

    std::size_t a = line.find_first_of(' ');
    std::size_t b = line.find_last_of(' ');

    if ((a == std::string::npos) || (a == b) || (b == line.size() - 1))
    {
        LOGTHROW(err4, std::runtime_error)
            << "Error parsing property line: " << line;
    }

    prop.type = line.substr(a + 1, b - a - 1);
    prop.name = line.substr(b + 1, line.size() - b - 1);
    return prop;
}

/// Parses lines of one element (vertex, face, ...) in header
/// After return, line contains a beginning of the next property or "end_header"
PlyElement parseElement(std::ifstream& f, std::string& line)
{
    std::string keyword, name;
    std::size_t count;

    if (!(std::istringstream(line) >> keyword >> name >> count))
    {
        LOGTHROW(err4, std::runtime_error)
            << "Error parsing element defined by: " << line;
    }

    if (keyword != "element")
    {
        LOGTHROW(err4, std::runtime_error)
            << "Error parsing element - does not start with \"element\": "
            << line;
    }

    PlyElement element;
    element.name = name;
    element.count = count;

    std::getline(f, line);
    ba::trim(line);
    throwOnEof(f);

    while (ba::starts_with(line, "property"))
    {
        element.props.push_back(parseProperty(line));

        std::getline(f, line);
        ba::trim(line);
        throwOnEof(f);
    }

    if (element.props.empty())
    {
        LOGTHROW(err4, std::runtime_error)
            << "Encountered element with no properties: " << name;
    }

    return element;
}

/// Parses the header section of PLY file and returns a vector of elements with their properties
std::vector<PlyElement> parseHeader(std::ifstream& f)
{
    std::string line;

    std::getline(f, line);
    ba::trim(line);
    throwOnEof(f);

    std::vector<PlyElement> elements;

    if (line != "ply")
    {
        LOGTHROW(err4, std::runtime_error)
            << "Unexpected start of PLY file: " << line;
    }

    do
    {
        if (ba::starts_with(line, "element"))
        {
            elements.push_back(parseElement(f, line));
        }
        else
        {
            std::getline(f, line);
            ba::trim(line);
            throwOnEof(f);
        }
    }
    while (!ba::starts_with(line, "end_header"));

    return elements;
}


void parsePly(PlyParserBase& parser, const fs::path& path)
{
    LOG(info1) << "Loading PLY file from " << path << ".";
    std::ifstream f(path.native());
    if (!f.good()) {
        LOGTHROW(err4, std::runtime_error)
            << "Unable to open file " << path << ">.";
    }
    // not setting failbit - is set on eof
    f.exceptions(std::ios::badbit);

    // load header
    std::vector<PlyElement> elements = parseHeader(f);
    parser.loadHeader(elements);
    
    // get face/vertex range
    std::size_t fStart = 0;
    std::size_t fEnd = 0;
    std::size_t vStart = 0;
    std::size_t vEnd = 0;
    
    std::size_t idx = 0;
    for (const auto& el : elements)
    {
        if (el.name == "face")
        {
            if (fEnd != 0)
            {
                LOGTHROW(err4, std::runtime_error)
                    << "Multiple elements named \"face\" found in header.";
            }
            fStart = idx;
            fEnd = idx + el.count;
        }
        else if (el.name == "vertex")
        {
            if (vEnd != 0)
            {
                LOGTHROW(err4, std::runtime_error)
                    << "Multiple elements named \"vertex\" found in header.";
            }
            vStart = idx;
            vEnd = idx + el.count;
        }
        idx += el.count;
    }
    std::size_t totNum = idx;  // expected number of lines

    LOG(info1) << "Expecting vertices on lines " << vStart << " - " << vEnd
               << "; faces on lines " << fStart << " - " << fEnd;

    // read the rest of the file
    idx = 0;
    std::size_t fRead = 0;
    std::size_t vRead = 0;
    while (!f.eof())
    {
        f.peek(); // skip last line at EOF
        if (f.eof()) { break; }

        if (idx >= fStart && idx < fEnd)
        {
            parser.addFace(f);
            ++fRead;
        }
        else if (idx >= vStart && idx < vEnd)
        {
            parser.addVertex(f);
            ++vRead;
        }
        else
        {
            parser.addOtherElement(f, idx);
        }

        ++idx;

        // deal with newline character
        char c;
        f.get(c);
        if (f.eof()) { break; } // no newline at the end of file
        if ((c != '\r') && (c != '\n'))
        {
            LOGTHROW(err4, std::runtime_error)
                << "Parser did not reach the EOL - next character: " << c;
        }
        if (c == '\r') // deal with "\r\n" or "\r" eol
        {
            c = f.peek();
            if (c == '\n') { f.get(c); }
        }
    }

    // check that the number of loaded elements corresponds to the header
    if ((fEnd - fStart) != fRead)
    {
        LOG(warn4)
            << "Expected " << fEnd - fStart << " faces (from header) but found "
            << fRead;
    }
    if ((vEnd - vStart) != vRead)
    {
        LOG(warn4)
            << "Expected " << vEnd - vStart
            << " vertices (from header) but found " << vRead;
    }
    if (idx != totNum)
    {
        LOG(warn4)
            << "Expected " << totNum << " lines (from header) but found "
            << idx;
    }

    f.close();
}

}  // namespace geometry
