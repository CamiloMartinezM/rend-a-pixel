#include <lightwave.hpp>

namespace lightwave {

    class Diffuse : public Bsdf {
        ref<Texture> m_albedo;

        public:

        Diffuse(const Properties& properties) {
            m_albedo = properties.get<Texture>("albedo");
        }

        BsdfEval evaluate(const Point2& uv, const Vector& wo, const Vector& wi) const override {
            NOT_IMPLEMENTED; // This is for assignment 3
        }

        BsdfSample sample(const Point2& uv, const Vector& wo,
                          Sampler& rng) const override {
            // Sample a random point
            const Point2 sampledPoint = Point2(rng.next(), rng.next());
            
            // Sample a direction on the cosine-weighted hemisphere
            Vector sampledVector = squareToCosineHemisphere(sampledPoint);

            // Evaluate the BSDF for the sampled direction
            Color bsdfValue = m_albedo->evaluate(uv);

            // Compute the probability of sampling the direction
            // float pdf = cosineHemispherePdf(sampledVector);

            // Correct for the foreshortening term cos 𝜔𝑖
            // if (pdf != 0)
            //    bsdfValue *= std::abs(Frame::cosTheta(sampledVector)) / pdf;

            return BsdfSample(sampledVector.normalized(), bsdfValue);
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