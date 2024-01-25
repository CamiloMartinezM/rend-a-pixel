#include "fresnel.hpp"
#include "microfacet.hpp"
#include <lightwave.hpp>

namespace lightwave
{
    class RoughConductor : public Bsdf
    {
      private:
        ref<Texture> m_reflectance;
        ref<Texture> m_roughness;

      public:
        RoughConductor(const Properties &properties)
        {
            m_reflectance = properties.get<Texture>("reflectance");
            m_roughness = properties.get<Texture>("roughness");
        }

        BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi) const override
        {
            // Using the squared roughness parameter results in a more gradual
            // transition from specular to rough. For numerical stability, we avoid
            // extremely specular distributions (alpha values below 10^-3)
            const auto alpha = std::max(float(1e-3), sqr(m_roughness->scalar(uv)));

            // hints:
            // * the microfacet normal can be computed from `wi' and `wo'
            Vector nWi = wi.normalized();
            Vector nWo = wo.normalized();
            Vector wm = (wi + wo).normalized();
            float D = microfacet::evaluateGGX(alpha, wm);
            float G1wi = microfacet::smithG1(alpha, wm, nWi);
            float G1wo = microfacet::smithG1(alpha, wm, nWo);
            Color R = m_reflectance->evaluate(uv);

            // Frame::absCosTheta(wi) cancels out from the denominator
            float denominator = 4 * Frame::absCosTheta(nWo);

            return {.value = R * D * G1wi * G1wo / denominator, .pdf = microfacet::pdfGGXVNDF(alpha, wm, nWo)};
        }

        BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng) const override
        {
            const auto alpha = std::max(float(1e-3), sqr(m_roughness->scalar(uv)));

            // hints:
            // * do not forget to cancel out as many terms from your equations as possible!
            //   (the resulting sample weight is only a product of two factors)
            Vector nWo = wo.normalized();
            Vector m = microfacet::sampleGGXVNDF(alpha, nWo, rng.next2D()); // Sample microfacet normal m
            Vector wi = reflect(wo, m).normalized(); // Reflect wo about m to get outgoing direction wi
            Vector wm = (wi + wo).normalized();
            float G1wi = microfacet::smithG1(alpha, wm, wi);
            Color weight = m_reflectance->evaluate(uv) * G1wi;
            return {.wi = wi, .weight = weight, .pdf = microfacet::pdfGGXVNDF(alpha, wm, nWo)};
        }

        std::string toString() const override
        {
            return tfm::format("RoughConductor[\n"
                               "  reflectance = %s,\n"
                               "  roughness = %s\n"
                               "]",
                               indent(m_reflectance), indent(m_roughness));
        }
    };
} // namespace lightwave

REGISTER_BSDF(RoughConductor, "roughconductor")