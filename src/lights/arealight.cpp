#include <lightwave.hpp>

namespace lightwave
{
    class AreaLight final : public Light
    {
      private:
        // Pointer to the Instance associated with the area light
        ref<Instance> m_shape;

        /**
         * @brief Computes the DirectLightSample based on the provided origin points and sampledArea of the shape.
         * @param origin Origin point.
         * @param sampledArea The result of sampling an AreaSample of the shape. Every shape has its way of doing this.
         * @return DirectLightSample instance.
         */
        DirectLightSample computeDirectLightSample(const Point &origin, const AreaSample &sampledArea) const
        {
            // Vector from the sampled point on the light source towards the origin
            // By convention, vectors point out of the surface
            Vector sampledAreaWo = origin - sampledArea.position;

            if (sampledArea.pdf == 0 || sampledAreaWo.lengthSquared() == 0)
                return DirectLightSample::invalid();

            Vector wi = (-sampledAreaWo).normalized();
            float distance = sampledAreaWo.length();

            // Vector on the local coordinate system
            Vector localSampledAreaWo = sampledArea.frame.toLocal(sampledAreaWo).normalized();

            // Evaluate emission at the sampled point on the shape's surface UV coordinates
            Color emission = m_shape->emission()->evaluate(sampledArea.uv, localSampledAreaWo).value;

            // Contribution of the angle between the surface normal pointing out and the incoming ray angle pointing
            // from the sampled point on the area light towards the origin point
            float cosTheta = Frame::absCosTheta(localSampledAreaWo);

            // Adjust intensity based on the area of the light (probability) and the angle of the surface normal and
            // the incoming ray angle with respect to the surface normal
            Color intensity = (emission * cosTheta) / (sampledArea.pdf * sqr(distance));

            return DirectLightSample{.wi = wi, .weight = intensity, .distance = distance};
        }

      public:
        AreaLight(const Properties &properties)
        {
            m_shape = properties.getChild<Instance>("instance");
            m_shape->setLight(this);
        }

        DirectLightSample sampleDirect(const Point &origin, Sampler &rng) const override
        {
            const AreaSample sampledArea = m_shape->sampleArea(rng);
            return computeDirectLightSample(origin, sampledArea);
        }

        DirectLightSample sampleDirect(const Point &origin, Sampler &rng, const SurfaceEvent &ref) const override
        {
            const AreaSample sampledArea = m_shape->sampleArea(rng, ref);
            return computeDirectLightSample(origin, sampledArea);
        }

        bool canBeIntersected() const override
        {
            return m_shape->isVisible();
        }

        std::string toString() const override
        {
            return tfm::format("AreaLight[\n"
                               "  shape = %s\n"
                               "]",
                               indent(m_shape));
        }
    };
} // namespace lightwave

REGISTER_LIGHT(AreaLight, "area")
