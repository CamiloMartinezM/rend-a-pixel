#include <lightwave.hpp>

namespace lightwave
{
    class Diffuse : public Bsdf
    {
        ref<Texture> m_albedo;

      public:
        Diffuse(const Properties &properties)
        {
            m_albedo = properties.get<Texture>("albedo");
        }

        BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi) const override
        {
            Vector wiCorrected = wi;

            // Check if wo and wi are in different hemispheres
            if (!Frame::sameHemisphere(wo, wi))
            {
                // Flip the incoming direction wi to the same hemisphere as wo
                wiCorrected = Vector(wi.x(), wi.y(), -wi.z());
            }

            return {m_albedo->evaluate(uv) * Frame::absCosTheta(wiCorrected) / Pi};
        }

        BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng) const override
        {
            Color albedo = m_albedo->evaluate(uv);
            Vector wi = squareToCosineHemisphere(rng.next2D());

            // Check if wo is in the opposite direction of the surface normal
            // If it is, flip wi to the opposite hemisphere
            wi *= Vector(1.0f, 1.0f, wo.z() < 0 ? -1.0f : 1.0f);
            return BsdfSample{wi, albedo};
        }

        std::string toString() const override
        {
            return tfm::format("Diffuse[\n"
                               "  albedo = %s\n"
                               "]",
                               indent(m_albedo));
        }
    };
} // namespace lightwave

REGISTER_BSDF(Diffuse, "diffuse")