#include <lightwave.hpp>

#include "../core/plyparser.hpp"
#include "accel.hpp"

namespace lightwave {

/**
 * @brief A shape consisting of many (potentially millions) of triangles, which share an index and vertex buffer.
 * Since individual triangles are rarely needed (and would pose an excessive amount of overhead), collections of
 * triangles are combined in a single shape.
 */
class TriangleMesh : public AccelerationStructure {
    /**
     * @brief The index buffer of the triangles.
     * The n-th element corresponds to the n-th triangle, and each component of the element corresponds to one
     * vertex index (into @c m_vertices ) of the triangle.
     * This list will always contain as many elements as there are triangles.
     */
    std::vector<Vector3i> m_triangles;
    /**
     * @brief The vertex buffer of the triangles, indexed by m_triangles.
     * Note that multiple triangles can share vertices, hence there can also be fewer than @code 3 * numTriangles @endcode
     * vertices.
     */
    std::vector<Vertex> m_vertices;
    /// @brief The file this mesh was loaded from, for logging and debugging purposes.
    std::filesystem::path m_originalPath;
    /// @brief Whether to interpolate the normals from m_vertices, or report the geometric normal instead.
    bool m_smoothNormals;

protected:
    int numberOfPrimitives() const override {
        return int(m_triangles.size());
    }

    bool intersect(int primitiveIndex, const Ray &ray, Intersection &its, Sampler &rng) const override {
     //   NOT_IMPLEMENTED

        // hints:
        // * use m_triangles[primitiveIndex] to get the vertex indices of the triangle that should be intersected
        // * if m_smoothNormals is true, interpolate the vertex normals from m_vertices
        //   * make sure that your shading frame stays orthonormal!
        // * if m_smoothNormals is false, use the geometrical normal (can be computed from the vertex positions)
    const Vector3i &triIndices = m_triangles[primitiveIndex];
    const Vertex &v0 = m_vertices[triIndices[0]];
    const Vertex &v1 = m_vertices[triIndices[1]];
    const Vertex &v2 = m_vertices[triIndices[2]];

    Vector edge1 = v1.position - v0.position;
    Vector edge2 = v2.position - v0.position;
    Vector pvec = ray.direction.cross(edge2);
    float det = edge1.dot(pvec);

    if (std::abs(det) < Epsilon) return false;

    float invDet = 1.0f / det;

    Vector tvec = ray.origin - v0.position;
    float u = tvec.dot(pvec) * invDet;
    if (u < 0 || u > 1) return false;

    Vector qvec = tvec.cross(edge1);
    float v = ray.direction.dot(qvec) * invDet;
    if (v < 0 || u + v > 1) return false;

    float t = edge2.dot(qvec) * invDet;

    if (t < 0) return false;

    // If m_smoothNormals is true, interpolate the normals
    if (m_smoothNormals) {
        //float w = 1.0f - u - v;
        Vector2 bary{u,v};
        its.frame.normal =Vertex::interpolate(bary, v0, v1, v2).normal.normalized();
    } else {
        // Compute plane's normal. If the dot product of normal and ray direction is positive,
        // the triangle is backfacing, which we don't want to hit
        its.frame.normal = edge1.cross(edge2).normalized();
       //////// if (its.frame.normal.dot(ray.direction) > 0.0f) its.frame.normal = -its.frame.normal; // Flip the normal if necessary
    }

    // Compute hit position and update intersection record
   // its.t = t;
   // its.p = ray(t); // Assuming ray(t) correctly computes the point on the ray at t
    // Assuming that the Vertex structure contains UV coordinates
   // its.uv = Vector(u, v); 

    return true;



    }

    Bounds getBoundingBox(int primitiveIndex) const override {
       // NOT_IMPLEMENTED
        const Vector3i &tri = m_triangles[primitiveIndex];
    const Vertex &v0 = m_vertices[tri[0]];
    const Vertex &v1 = m_vertices[tri[1]];
    const Vertex &v2 = m_vertices[tri[2]];

    float minX = std::min({v0.position.x(), v1.position.x(), v2.position.x()});
    float minY = std::min({v0.position.y(), v1.position.y(), v2.position.y()});
    float minZ = std::min({v0.position.z(), v1.position.z(), v2.position.z()});

    float maxX = std::max({v0.position.x(), v1.position.x(), v2.position.x()});
    float maxY = std::max({v0.position.y(), v1.position.y(), v2.position.y()});
    float maxZ = std::max({v0.position.z(), v1.position.z(), v2.position.z()});

    Point minPoint(minX, minY, minZ);
    Point maxPoint(maxX, maxY, maxZ);

    return Bounds(minPoint, maxPoint);
    }

    Point getCentroid(int primitiveIndex) const override {
        //NOT_IMPLEMENTED
        const Vector3i &tri = m_triangles[primitiveIndex];
    const Vertex &v0 = m_vertices[tri[0]];
    const Vertex &v1 = m_vertices[tri[1]];
    const Vertex &v2 = m_vertices[tri[2]];

    // Sum up each component of the vertices
    float sumX = v0.position.x() + v1.position.x() + v2.position.x();
    float sumY = v0.position.y() + v1.position.y() + v2.position.y();
    float sumZ = v0.position.z() + v1.position.z() + v2.position.z();

    // Compute the average for each component
    float centerX = sumX / 3.0f;
    float centerY = sumY / 3.0f;
    float centerZ = sumZ / 3.0f;

    return Point(centerX, centerY, centerZ);
    }

public:
    TriangleMesh(const Properties &properties) {
        m_originalPath = properties.get<std::filesystem::path>("filename");
        m_smoothNormals = properties.get<bool>("smooth", true);
        readPLY(m_originalPath.string(), m_triangles, m_vertices);
        logger(EInfo, "loaded ply with %d triangles, %d vertices",
            m_triangles.size(),
            m_vertices.size()
        );
        buildAccelerationStructure();
    }

    AreaSample sampleArea(Sampler &rng) const override {
        // only implement this if you need triangle mesh area light sampling for your rendering competition
        NOT_IMPLEMENTED
    }

    std::string toString() const override {
        return tfm::format(
            "Mesh[\n"
            "  vertices = %d,\n"
            "  triangles = %d,\n"
            "  filename = \"%s\"\n"
            "]",
            m_vertices.size(),
            m_triangles.size(),
            m_originalPath.generic_string()
        );
    }
};

}

REGISTER_SHAPE(TriangleMesh, "mesh")
