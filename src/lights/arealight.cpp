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
        DirectLightSample computeDirectLightSample(const Point &origin, AreaSample sampledArea) const
        {
            Vector sampledAreaVector = sampledArea.position - origin; // Vector from the origin to the sampled point

            if (sampledArea.pdf == 0 || sampledAreaVector.lengthSquared() == 0)
                return DirectLightSample::invalid();

            Vector direction = sampledAreaVector.normalized();
            float distance = sampledAreaVector.length();

            // Evaluate emission at the sampled point on the shape's surface UV coordinates
            Color emission = m_shape->emission()->evaluate(sampledArea.uv, sampledArea.frame.toLocal(direction)).value;

            // Adjust intensity based on the area of the light (probability), not diving by squared distance, therefore
            // using radiance as a measure
            Color intensity = emission / sampledArea.pdf;

            return DirectLightSample{.wi = direction, .weight = intensity, .distance = distance};
        }

      public:
        AreaLight(const Properties &properties)
        {
            m_shape = properties.getChild<Instance>("instance");
        }

        DirectLightSample sampleDirect(const Point &origin, Sampler &rng) const override
        {
            AreaSample sampledArea = m_shape->sampleArea(rng);
            return computeDirectLightSample(origin, sampledArea);
        }

        DirectLightSample sampleDirect(const Point &origin, Sampler &rng, const Intersection &ref) const override
        {
            AreaSample sampledArea = m_shape->sampleArea(rng, ref);
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
