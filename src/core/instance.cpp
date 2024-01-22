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

        // If m_flipNormal is true, flip the direction of the bitangent
        transformedBitangent *= m_flipNormal ? -1 : 1;

        // Normalize the tangent and bitangent to ensure the frame remains orthonormal
        surf.frame.tangent = transformedTangent.normalized();
        surf.frame.bitangent = transformedBitangent.normalized();
        surf.frame.normal = surf.frame.tangent.cross(surf.frame.bitangent).normalized();

        // Adjust PDF for transformed area, measured by the cross product, to account for the transforms which can
        // affect the sampling density
        surf.pdf /= abs(crossProduct.length());
    }

    void Instance::inverseTransformFrame(SurfaceEvent &surf) const
    {
        surf.position = m_transform->inverse(surf.position);
        Vector transformedTangent = m_transform->inverse(surf.frame.tangent);
        Vector transformedBitangent = m_transform->inverse(surf.frame.bitangent);

        // The length of the cross product gives us the area spanned by tangent and bitangent, which is a measure of the
        // change in surface area, which can affect the sampling density
        Vector crossProduct = transformedTangent.cross(transformedBitangent);

        // If m_flipNormal is true, flip the direction of the bitangent
        transformedBitangent *= m_flipNormal ? -1 : 1;

        // Normalize the tangent and bitangent to ensure the frame remains orthonormal
        surf.frame.tangent = transformedTangent.normalized();
        surf.frame.bitangent = transformedBitangent.normalized();
        surf.frame.normal = surf.frame.tangent.cross(surf.frame.bitangent).normalized();

        // Adjust PDF for transformed area, measured by the cross product, to account for the transforms which can
        // affect the sampling density
        surf.pdf /= abs(crossProduct.length());
    }

    bool Instance::intersect(const Ray &worldRay, Intersection &its, Sampler &rng) const
    {
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
        inverseTransformFrame(localRef);
        AreaSample sample = m_shape->sampleArea(rng, localRef);
        transformFrame(sample);
        return sample;
    }

    inline float Instance::sampledDirectionPdf(const Vector &sampledVector) const
    {
        return m_shape->sampledDirectionPdf(sampledVector);
    }

} // namespace lightwave

REGISTER_CLASS(Instance, "instance", "default")
