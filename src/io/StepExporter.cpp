// ------------------------------------------------------------------------------
// Project: BladeGen
// Copyright(c) 2026, Onur Tuncer, PhD, Istanbul Technical University
//
// SPDX-License-Identifier: GPL-3.0-or-later
// License-Filename: LICENSE
// ------------------------------------------------------------------------------
/// @file io/StepExporter.cpp
/// @brief STEP export implementation.

#include "io/StepExporter.hpp"

#include <stdexcept>
#include <string>

// OCCT
#include <TopoDS_Shape.hxx>
#include <STEPControl_Writer.hxx>
#include <STEPControl_StepModelType.hxx>
#include <Interface_Static.hxx>
#include <IFSelect_ReturnStatus.hxx>

namespace BladeGen::IO {

StepExporter::StepExporter(StepProtocol protocol)
    : m_protocol(protocol)
{}

void StepExporter::AddShape(const TopoDS_Shape& shape, const std::string& label)
{
    m_shapes.emplace_back(shape, label);
}

void StepExporter::Write(const std::filesystem::path& path) const
{
    if (m_shapes.empty())
        throw std::runtime_error("StepExporter::Write — no shapes staged");

    STEPControl_Writer writer;

    // Set the STEP schema
    const char* schema = (m_protocol == StepProtocol::AP203) ? "AP203" : "AP214IS";
    Interface_Static::SetCVal("write.step.schema", schema);

    for (const auto& [shape, label] : m_shapes) {
        if (!label.empty())
            Interface_Static::SetCVal("write.step.product.name", label.c_str());

        const IFSelect_ReturnStatus status =
            writer.Transfer(shape, STEPControl_AsIs);

        if (status != IFSelect_RetDone)
            throw std::runtime_error("StepExporter: OCCT Transfer failed for shape '" + label + "'");
    }

    const IFSelect_ReturnStatus writeStatus =
        writer.Write(path.string().c_str());

    if (writeStatus != IFSelect_RetDone)
        throw std::runtime_error("StepExporter: OCCT Write failed for '" + path.string() + "'");
}

} // namespace BladeGen::IO
