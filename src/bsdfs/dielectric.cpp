#include "fresnel.hpp"
#include <lightwave.hpp>

namespace lightwave {

    class Dielectric : public Bsdf {
        ref<Texture> m_ior;
        ref<Texture> m_reflectance;
        ref<Texture> m_transmittance;

        public:
        Dielectric(const Properties& properties) {
            m_ior = properties.get<Texture>("ior");
            m_reflectance = properties.get<Texture>("reflectance");
            m_transmittance = properties.get<Texture>("transmittance");
        }

        BsdfEval evaluate(const Point2& uv, const Vector& wo,
                          const Vector& wi) const override {
            // the probability of a light sample picking exactly the direction `wi'
            // that results from reflecting or refracting `wo' is zero, hence we can
            // just ignore that case and always return black
            return BsdfEval::invalid();
        }

        BsdfSample sample(const Point2& uv, const Vector& wo,
                          Sampler& rng) const override {
            float ior = m_ior->scalar(uv);
            bool entering = Frame::cosTheta(wo) > 0;
            Vector normal = Vector(0.f, 0.f, 1.f) * (entering ? 1.f : -1.f); // Flip the normal for exiting rays
            float eta = entering ? ior : 1.f / ior;
            float F = fresnelDielectric(entering ? Frame::cosTheta(wo) : -Frame::cosTheta(wo), eta);

            if (rng.next() < F) {
                Vector wi = reflect(wo, normal);
                Color reflectance = m_reflectance->evaluate(uv);
                return BsdfSample(wi, reflectance * F);
            }
            else {
                Vector wi = refract(wo, normal, eta);
                if (!wi.isZero()) {
                    float cosThetaT = Frame::absCosTheta(wi);
                    Color transmittance = m_transmittance->evaluate(uv) * cosThetaT * 1 / (eta * eta);
                    return BsdfSample(wi, transmittance * (1.f - F));
                }

                return BsdfSample(reflect(wo, normal), m_reflectance->evaluate(uv)); // Handle total internal reflection
            }
        }

        std::string toString() const override {
            return tfm::format("Dielectric[\n"
                               "  ior           = %s,\n"
                               "  reflectance   = %s,\n"
                               "  transmittance = %s\n"
                               "]",
                               indent(m_ior), indent(m_reflectance),
                               indent(m_transmittance));
        }
    };

} // namespace lightwave

REGISTER_BSDF(Dielectric, "dielectric")
