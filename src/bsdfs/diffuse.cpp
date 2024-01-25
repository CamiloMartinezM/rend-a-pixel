#include "diffuse.hpp"
#include <lightwave.hpp>

namespace lightwave
{
    class Diffuse : public Bsdf
    {
      private:
        ref<Texture> m_albedo;

      public:
        Diffuse(const Properties &properties)
        {
            m_albedo = properties.get<Texture>("albedo");
        }

        BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi) const override
        {
            Vector nWi = wi.normalized();
            float cosTheta = Frame::sameHemisphere(wo, wi) ? Frame::absCosTheta(nWi) : max(0.0f, Frame::cosTheta(nWi));
            return {.value = m_albedo->evaluate(uv) * cosTheta * InvPi, .pdf = diffuse::pdf(wo, nWi)};
        }

        BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng) const override
        {
            Color albedo = m_albedo->evaluate(uv);
            Vector wi = squareToCosineHemisphere(rng.next2D());

            // Check if wo is in the opposite direction of the surface normal. If it is, flip wi to the opposite 
            // hemisphere, i.e, to make it the same as wo
            if (!Frame::sameHemisphere(wo, wi))
                wi *= Vector(1.0f, 1.0f, -1.0f);
            
            return {.wi = wi, .weight = albedo, .pdf = diffuse::pdf(wo, wi)};
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