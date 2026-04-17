/// @file io/IgesExporter.cpp
/// @brief IGES export implementation.

#include "io/IgesExporter.hpp"

#include <stdexcept>

#include <IGESControl_Writer.hxx>
#include <IFSelect_ReturnStatus.hxx>
#include <Interface_Static.hxx>
#include <TopoDS_Shape.hxx>

namespace PCAD::IO {

void IgesExporter::AddShape(const TopoDS_Shape& shape, const std::string& label)
{
    m_shapes.emplace_back(shape, label);
}

void IgesExporter::Write(const std::filesystem::path& path) const
{
    if (m_shapes.empty())
        throw std::runtime_error("IgesExporter::Write — no shapes staged");

    // BRep mode (0): preserves OCCT topology faithfully.
    // Analytic mode (1): converts to IGES analytic entities — less portable.
    IGESControl_Writer writer("MM", 0);

    for (const auto& [shape, label] : m_shapes) {
        if (!label.empty())
            Interface_Static::SetCVal("write.iges.header.product", label.c_str());

        if (!writer.AddShape(shape))
            throw std::runtime_error(
                "IgesExporter: OCCT AddShape failed for '" + label + "'");
    }

    writer.ComputeModel();

    if (!writer.Write(path.string().c_str()))
        throw std::runtime_error(
            "IgesExporter: OCCT Write failed for '" + path.string() + "'");
}

} // namespace PCAD::IO
