#include <lightwave.hpp>

namespace lightwave
{
    class Conductor : public Bsdf
    {
      private:
        ref<Texture> m_reflectance;

      public:
        Conductor(const Properties &properties)
        {
            m_reflectance = properties.get<Texture>("reflectance");
        }

        BsdfEval evaluate(const Point2 &uv, const Vector &wo, const Vector &wi) const override
        {
            // The probability of a light sample picking exactly the direction `wi' that results from reflecting or
            // refracting `wo' is zero, hence we can just ignore that case and always return black and pdf = 0.0f
            return BsdfEval::invalid();
        }

        BsdfSample sample(const Point2 &uv, const Vector &wo, Sampler &rng) const override
        {
            Vector wi = reflect(wo, Vector(0.0f, 0.0f, 1.0f)).normalized(); // Reflect wo about the normal to get wi
            Color reflectanceValue = m_reflectance->evaluate(uv); // Compute the BSDF reflectance value

            // There is no randomness in the conductor reflection, so the pdf is a delta function and
            // the weight is the reflectance value itself
            return {.wi = wi, .weight = reflectanceValue, .pdf = Infinity};
        }

        std::string toString() const override
        {
            return tfm::format("Conductor[\n"
                               "  reflectance = %s\n"
                               "]",
                               indent(m_reflectance));
        }
    };

} // namespace lightwave

REGISTER_BSDF(Conductor, "conductor")
