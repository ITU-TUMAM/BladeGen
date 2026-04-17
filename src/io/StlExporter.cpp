/// @file io/StlExporter.cpp
/// @brief STL export implementation.

#include "io/StlExporter.hpp"

#include <stdexcept>
#include <string>

// OCCT
#include <TopoDS_Shape.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <StlAPI_Writer.hxx>

namespace PCAD::IO {

void StlExporter::Write(const TopoDS_Shape& shape,
                        const std::filesystem::path& path) const
{
    if (shape.IsNull())
        throw std::runtime_error("StlExporter::Write — shape is null");

    // Mesh the shape first — StlAPI_Writer produces an empty file on unmeshed shapes.
    BRepMesh_IncrementalMesh mesher(shape, m_linearDeflection,
                                    /*isRelative=*/false,
                                    m_angularDeflection);
    mesher.Perform();

    if (!mesher.IsDone())
        throw std::runtime_error("StlExporter: BRepMesh failed for '" + path.string() + "'");

    StlAPI_Writer writer;
    const bool ok = writer.Write(shape, path.string().c_str());

    if (!ok)
        throw std::runtime_error("StlExporter: StlAPI_Writer failed for '" + path.string() + "'");
}

} // namespace PCAD::IO
