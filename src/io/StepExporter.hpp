/// @file io/StepExporter.hpp
/// @brief STEP file export via OpenCASCADE STEPControl_Writer.
///
/// Stages one or more shapes, then writes a single STEP file.
/// Supports AP203 (simple geometry) and AP214 (colours + layers).

#pragma once

#include <filesystem>
#include <string>
#include <vector>

class TopoDS_Shape;

namespace PCAD::IO {

enum class StepProtocol { AP203, AP214 };

class StepExporter
{
public:
    explicit StepExporter(StepProtocol protocol = StepProtocol::AP214);

    /// Stage a shape with an optional label (used in the STEP product name).
    void AddShape(const TopoDS_Shape& shape, const std::string& label = "");

    /// Write all staged shapes to @p path.
    /// @throws std::runtime_error on OCCT write failure.
    void Write(const std::filesystem::path& path) const;

    void Clear() noexcept { m_shapes.clear(); }

private:
    StepProtocol                                   m_protocol;
    std::vector<std::pair<TopoDS_Shape, std::string>> m_shapes;
};

} // namespace PCAD::IO
