/// @file io/IgesExporter.hpp
/// @brief IGES file export via OpenCASCADE IGESControl_Writer.
///
/// Stages one or more shapes, then writes a single IGES file.
/// IGES is widely supported by legacy CAD and CFD mesh tools.

#pragma once

#include <filesystem>
#include <string>
#include <vector>

class TopoDS_Shape;

namespace PCAD::IO {

class IgesExporter
{
public:
    /// Stage a shape with an optional label (written to the IGES entity name).
    void AddShape(const TopoDS_Shape& shape, const std::string& label = "");

    /// Write all staged shapes to @p path.
    /// @throws std::runtime_error on OCCT write failure.
    void Write(const std::filesystem::path& path) const;

    void Clear() noexcept { m_shapes.clear(); }

private:
    std::vector<std::pair<TopoDS_Shape, std::string>> m_shapes;
};

} // namespace PCAD::IO
