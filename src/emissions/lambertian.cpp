#include <lightwave.hpp>

namespace lightwave {

    class Lambertian : public Emission {
        ref<Texture> m_emission;

        public:
        Lambertian(const Properties& properties) {
            m_emission = properties.get<Texture>("emission");
        }

        EmissionEval evaluate(const Point2& uv, const Vector& wo) const override {
            // Evaluate the texture at the given UV coordinates to get the emission color
            // - Check the sign of the cosine to not emit any light in the backwards direction of the lightsource. 
            // - Check whether the sign of your cosine in the Lambertian emission evaluation is negative and 
            //   return zero in that case.
            Color emissionColor = (Frame::cosTheta(wo) < 0)? Color(0.0f) : m_emission->evaluate(uv);
            
            // No π factor needed for the value.
            // The PDF for a uniformly emitting Lambertian surface
            return EmissionEval{ .value = emissionColor, .pdf = Inv2Pi }; 
        }

        std::string toString() const override {
            return tfm::format("Lambertian[\n"
                               "  emission = %s\n"
                               "]",
                               indent(m_emission));
        }
    };

} // namespace lightwave

REGISTER_EMISSION(Lambertian, "lambertian")
