#include <lightwave.hpp>
#include <math.h>

namespace lightwave {

    class Diffuse : public Bsdf {
        ref<Texture> m_albedo;

        public:

        Diffuse(const Properties& properties) {
            m_albedo = properties.get<Texture>("albedo");
        }

        BsdfEval evaluate(const Point2& uv, const Vector& wo, const Vector& wi) const override {
            Color albedo = m_albedo->evaluate(uv);
            Vector normalizedWo = wo.normalized();
            Vector normalizedWi = wi.normalized();

            // Compute the cosine of the angle between the normal and the incoming direction
            float cosTheta = normalizedWi.dot(normalizedWo.dot(Vector(0, 0, 1)) * Vector(0, 0, 1)); 
            Color brdfValue = albedo / M_PI * cosTheta;
            return { brdfValue };
        }

        BsdfSample sample(const Point2& uv, const Vector& wo,
                          Sampler& rng) const override {
            // Sample a direction on the cosine-weighted hemisphere
            Vector wi = squareToCosineHemisphere(uv);

            // Evaluate the BSDF for the sampled direction
            Color bsdfValue = m_albedo->evaluate(uv) / M_PI;

            // Compute the probability of sampling the direction
            float pdf = cosineHemispherePdf(wi);

            // Correct for the foreshortening term cos 𝜔𝑖
            bsdfValue *= std::abs(wi.y()) / pdf;

            // Return the sampled direction and the corresponding BSDF value
            return BsdfSample{ wi, bsdfValue };
        }

        std::string toString() const override {
            return tfm::format("Diffuse[\n"
                               "  albedo = %s\n"
                               "]",
                               indent(m_albedo));
        }
    };
}

REGISTER_BSDF(Diffuse, "diffuse")