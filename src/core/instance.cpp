#include <lightwave/core.hpp>
#include <lightwave/instance.hpp>
#include <lightwave/registry.hpp>
#include <lightwave/sampler.hpp>

namespace lightwave
{
    void Instance::transformFrame(SurfaceEvent &surf) const
    {
        // hints:
        // * transform the hitpoint and frame here
        // * if m_flipNormal is true, flip the direction of the bitangent (which in effect flips the normal)
        // * make sure that the frame is orthonormal (you are free to change the bitangent for this, but keep
        //   the direction of the transformed tangent the same)
        surf.position = m_transform->apply(surf.position);
        Vector transformedTangent = m_transform->apply(surf.frame.tangent);
        Vector transformedBitangent = m_transform->apply(surf.frame.bitangent);

        // The length of the cross product gives us the area spanned by tangent and bitangent, which is a measure of the
        // change in surface area, which can affect the sampling density
        Vector crossProduct = transformedTangent.cross(transformedBitangent);

        // Define a new world normal out of which a new frame will be built
        Vector worldNormal;
        if (m_normal) // Apply Normal Mapping feature, if it was provided for the instance
        {
            // Evaluate the uv coordinates to get the Color and convert it to a normal Vector
            Color localColor = m_normal->evaluate(surf.uv);
            Vector localNormal = Vector(localColor.r(), localColor.g(), localColor.b());

            // Map into [-1, 1]
            localNormal = 2.0f * localNormal - Vector(1.0f);
            worldNormal = m_transform->applyNormal(surf.frame.toWorld(localNormal));
        }
        else
        {
            // If m_flipNormal is true, flip the direction of the bitangent
            transformedBitangent *= m_flipNormal ? -1 : 1;
            worldNormal = transformedTangent.normalized().cross(transformedBitangent.normalized());
        }

        // Build an orthonormal frame with the new normal
        surf.frame = Frame(worldNormal.normalized());

        // Adjust PDF for transformed area, measured by the cross product, to account for the transforms which can
        // affect the sampling density
        surf.pdf /= abs(crossProduct.length());
    }

    void Instance::inverseTransformFrame(SurfaceEvent &surf) const
    {
        surf.position = m_transform->inverse(surf.position);
        Vector transformedTangent = m_transform->inverse(surf.frame.tangent);
        Vector transformedBitangent = m_transform->inverse(surf.frame.bitangent);

        // If m_flipNormal is true, flip the direction of the bitangent
        transformedBitangent *= m_flipNormal ? -1 : 1;

        // Normalize the tangent and bitangent to ensure the frame remains orthonormal
        Vector newNormal = transformedTangent.normalized().cross(transformedBitangent.normalized());
        surf.frame = Frame(newNormal.normalized());
    }

    bool Instance::intersect(const Ray &worldRay, Intersection &its, Sampler &rng) const
    {
        // Give the alpha mask property to the intersection object, so that the primitives can use it in their
        // intersect methods, i.e, Sphere::intersect, TriangleMesh::intersect, etc.
        if (m_alpha)
            its.alphaMask = m_alpha.get();

        if (!m_transform)
        {
            // fast path, if no transform is needed
            Ray localRay = worldRay;
            if (m_shape->intersect(localRay, its, rng))
            {
                its.instance = this;
                return true;
            }
            else
            {
                return false;
            }
        }

        const float previousT = its.t;
        Ray localRay;

        // hints:
        // * transform the ray (do not forget to normalize!)
        // * how does its.t need to change?
        localRay = m_transform->inverse(worldRay);
        const float localRayLength = localRay.direction.length();
        localRay = localRay.normalized();

        its.t = previousT * localRayLength;
        const bool wasIntersected = m_shape->intersect(localRay, its, rng);
        if (wasIntersected)
        {
            // hint: how does its.t need to change?
            its.t = its.t / localRayLength;
            its.instance = this;
            transformFrame(its);
            return true;
        }
        else
        {
            its.t = previousT;
            return false;
        }
    }

    Bounds Instance::getBoundingBox() const
    {
        if (!m_transform)
        {
            // fast path
            return m_shape->getBoundingBox();
        }

        const Bounds untransformedAABB = m_shape->getBoundingBox();
        if (untransformedAABB.isUnbounded())
        {
            return Bounds::full();
        }

        Bounds result;
        for (int point = 0; point < 8; point++)
        {
            Point p = untransformedAABB.min();
            for (int dim = 0; dim < p.Dimension; dim++)
            {
                if ((point >> dim) & 1)
                {
                    p[dim] = untransformedAABB.max()[dim];
                }
            }
            p = m_transform->apply(p);
            result.extend(p);
        }
        return result;
    }

    Point Instance::getCentroid() const
    {
        if (!m_transform)
        {
            return m_shape->getCentroid(); // fast path
        }
        return m_transform->apply(m_shape->getCentroid());
    }

    AreaSample Instance::sampleArea(Sampler &rng) const
    {
        AreaSample sample = m_shape->sampleArea(rng);
        transformFrame(sample);
        return sample;
    }

    AreaSample Instance::sampleArea(Sampler &rng, const SurfaceEvent &ref) const
    {
        SurfaceEvent localRef = ref;
        // localRef.position = m_transform->inverse(localRef.position);
        inverseTransformFrame(localRef);
        AreaSample sample = m_shape->sampleArea(rng, localRef);
        transformFrame(sample);
        return sample;
    }
} // namespace lightwave

REGISTER_CLASS(Instance, "instance", "default")
