#include <lightwave.hpp>

namespace lightwave {

    class Conductor : public Bsdf {
        ref<Texture> m_reflectance;

        public:
        Conductor(const Properties& properties) {
            m_reflectance = properties.get<Texture>("reflectance");
        }

        BsdfEval evaluate(const Point2& uv, const Vector& wo,
                          const Vector& wi) const override {
            // the probability of a light sample picking exactly the direction `wi'
            // that results from reflecting `wo' is zero, hence we can just ignore
            // that case and always return black
            return BsdfEval::invalid();
        }

        BsdfSample sample(const Point2& uv, const Vector& wo,
                          Sampler& rng) const override {
            Vector wi = reflect(wo, Vector(0.0f, 0.0f, 1.0f)); // Reflect wo about the normal to get wi
            Color reflectanceValue = m_reflectance->evaluate(uv); // Compute the BSDF reflectance value

            // There is no randomness in the conductor reflection, so the pdf is a delta function and 
            // the weight is the reflectance value itself
            return BsdfSample(wi, reflectanceValue);
        }

        std::string toString() const override {
            return tfm::format("Conductor[\n"
                               "  reflectance = %s\n"
                               "]",
                               indent(m_reflectance));
        }
    };

} // namespace lightwave

REGISTER_BSDF(Conductor, "conductor")
