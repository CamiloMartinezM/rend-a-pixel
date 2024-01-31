#include <lightwave.hpp>

namespace lightwave
{
    class AlbedoIntegrator : public SamplingIntegrator
    {
      public:
        AlbedoIntegrator(const Properties &properties) : SamplingIntegrator(properties)
        {
        }

        Color Li(const Ray &ray, Sampler &rng) override
        {
            Intersection its = m_scene->intersect(ray, rng);
            
            // Make sure that there is an intersection and it's a BSDF the hit surface
            if (!its || !its.instance->bsdf())
                return Color::black();

            const auto &bsdf = its.instance->bsdf();
            return bsdf->getAlbedo(its.uv);
        }

        std::string toString() const override
        {
            return "AlbedoIntegrator[]";
        }
    };
} // namespace lightwave

REGISTER_INTEGRATOR(AlbedoIntegrator, "albedo")
