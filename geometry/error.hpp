#ifndef geometry_error_hpp_included_
#define geometry_error_hpp_included_

#include <stdexcept>
#include <string>

namespace geometry {

#define DECLARE_EXCEPTION(type, base) \
    struct type : public base { type(const std::string &msg) : base(msg) {} }

DECLARE_EXCEPTION(Error, std::runtime_error);
DECLARE_EXCEPTION(BadFileFormat, Error);
DECLARE_EXCEPTION(VersionError, Error);

#undef DECLARE_EXCEPTION

} // namespace geometry

#endif // geometry_error_hpp_included_
