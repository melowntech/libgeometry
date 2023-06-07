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

#include <boost/algorithm/string.hpp>

namespace geometry
{

namespace ba = boost::algorithm;

/// Header parsing functions

/// Helper retrieving the next line of file, throws if EOF reached
inline void nextLine(std::ifstream& f, std::string& line)
{
    try
    {
        std::getline(f, line);
    }
    catch (const std::exception& e)
    {
        if (f.eof())
        {
            LOGTHROW(err4, std::runtime_error)
                << "Unexpected EOF reached while parsing PLY header";
        }
        throw e;
    }
    ba::trim(line);
}

/// Splits property line to name and type
bool parseProperty(const std::string& line, std::vector<PlyProperty>& props)
{
    if (!ba::starts_with(line, "property")) { return false; }

    std::size_t a = line.find_first_of(' ');
    std::size_t b = line.find_last_of(' ');
    if ((a == std::string::npos) || (a == b) || (b == line.size() - 1))
    {
        LOGTHROW(err4, std::runtime_error)
            << "Error parsing property line: " << line;
    }

    props.emplace_back(line.substr(a + 1, b - a - 1),          // prop type
                       line.substr(b + 1, line.size() - b - 1) // prop name

    );
    return true;
}

/// Checks if a line contains a comment and optionally parses it
bool parseComment(const std::string& line, std::vector<std::string>& comments)
{
    if (!ba::starts_with(line, "comment")) { return false; }

    comments.push_back(line.substr(7));
    ba::trim(comments.back()); // trim leading spaces
    return true;
}

/// Parses lines of one element (vertex, face, ...) in header
/// NB: After return, line contains a beginning of the next property, comment,
/// or "end_header"
bool parseElement(std::ifstream& f,
                  std::string& line,
                  std::vector<PlyElement>& elements,
                  std::vector<std::string>& comments)
{
    if (!ba::starts_with(line, "element")) { return false; }

    std::vector<std::string> tokens;
    ba::split(tokens, line, ba::is_space());
    if (tokens.size() != 3)
    {
        LOGTHROW(err4, std::runtime_error)
            << "Error parsing element definition line - expected three tokens: "
            << line;
    }

    const std::string& name = tokens[1];

    std::size_t count;
    try
    {
        count = std::stoul(tokens[2]);
    }
    catch (const std::invalid_argument& e)
    {
        LOGTHROW(err4, std::runtime_error)
            << "Error parsing element count from: " << tokens[2];
        throw;
    }
    catch (const std::out_of_range& e)
    {
        LOGTHROW(err4, std::runtime_error)
            << "Error parsing element count - element count out of range: "
            << tokens[2];
        throw;
    }

    elements.emplace_back(name, count);
    PlyElement& el = elements.back();
    // parse the element properties
    do
    {
        nextLine(f, line);
    }
    while (parseProperty(line, el.props) || parseComment(line, comments));

    if (el.props.empty())
    {
        LOGTHROW(err4, std::runtime_error)
            << "Encountered element with no properties: " << name;
    }

    return true;
}

/// Check "format" line in PLY file, throws on unexpected inputs
void checkFormatLine(std::ifstream& f)
{
    std::string line;
    nextLine(f, line); // read whole line

    std::vector<std::string> tokens;
    ba::split(tokens, line, ba::is_space());

    if (tokens.size() != 3)
    {
        LOGTHROW(err4, std::runtime_error)
            << "Error parsing format line - expected three tokens: "
            << line;
    }
    if (tokens[0] != "format")
    {
        LOGTHROW(err4, std::runtime_error)
            << "Expected the second line to contain \"format\" specification";
    }
    if (tokens[1] != "ascii")
    {
        LOGTHROW(err4, std::runtime_error)
            << "Only ASCII format is supported, got: " << tokens[1];
    }
    if (tokens[2] != "1.0")
    {
        LOGTHROW(err4, std::runtime_error)
            << "Only PLY version 1.0 supported, got: " << tokens[2];
    }
}

/// Parses the header section of PLY file and returns a vector of elements with
/// their properties and comments
std::pair<std::vector<PlyElement>, std::vector<std::string>>
    parseHeader(std::ifstream& f)
{
    // check magic number
    std::string line;
    nextLine(f, line);
    if (line != "ply")
    {
        LOGTHROW(err4, std::runtime_error)
            << "Unexpected magic number of PLY file: " << line;
    }

    // check line specifying PLY format ("format ascii 1.0")
    checkFormatLine(f);

    // parse elements & comments
    std::vector<PlyElement> elements;
    std::vector<std::string> comments;
    nextLine(f, line);
    do
    {
        if (parseElement(f, line, elements, comments))
        {
            continue; // parseElement() jumps to the line after element
                      // definition
        }
        else if (!parseComment(line, comments))
        {
            LOGTHROW(err4, std::runtime_error)
                << "Encountered unexpected line in header (not a comment nor "
                   "element definition): "
                << line;
        }

        nextLine(f, line);
    }
    while (!ba::starts_with(line, "end_header"));

    return std::make_pair(std::move(elements), std::move(comments));
}

void parsePly(PlyParserBase& parser, const fs::path& path)
{
    LOG(info1) << "Loading PLY file from " << path << ".";
    std::ifstream f(path.native());
    if (!f.good())
    {
        LOGTHROW(err4, std::runtime_error)
            << "Unable to open file " << path << ".";
    }
    f.exceptions(std::ios::badbit | std::ios::failbit);

    // load header
    std::vector<PlyElement> elements;
    std::vector<std::string> comments;
    std::tie(elements, comments) = parseHeader(f);
    parser.loadHeader(elements, comments);

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
    std::size_t totNum = idx; // expected number of lines

    LOG(info1) << "Expecting vertices on lines " << vStart << " - " << vEnd
               << "; faces on lines " << fStart << " - " << fEnd;

    // read the rest of the file
    for (idx = 0; idx < totNum; ++idx)
    {
        if (idx >= fStart && idx < fEnd)
        {
            parser.addFace(f);
        }
        else if (idx >= vStart && idx < vEnd)
        {
            parser.addVertex(f);
        }
        else
        {
            parser.addOtherElement(f, idx);
        }
    }

    f.close();
}

}  // namespace geometry
