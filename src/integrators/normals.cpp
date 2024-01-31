#include <lightwave.hpp>

namespace lightwave
{
    class NormalsIntegrator : public SamplingIntegrator
    {
      private:
        bool remap;

      public:
        NormalsIntegrator(const Properties &properties) : SamplingIntegrator(properties)
        {
            remap = properties.get<bool>("remap", true);
        }

        Color Li(const Ray &ray, Sampler &rng) override
        {
            Intersection its = m_scene->intersect(ray, rng);
            Vector normal = its.frame.normal;
            if (remap)
            {
                // Remapping between [-1, 1] to [0, 1]
                normal = (normal + Vector(1.0f)) * 0.5f;
            }
            return Color(normal);
        }

        std::string toString() const override
        {
            return tfm::format("CameraIntegrator[\n"
                               "  sampler = %s,\n"
                               "  image = %s,\n"
                               "]",
                               indent(m_sampler), indent(m_image));
        }
    };

} // namespace lightwave

REGISTER_INTEGRATOR(NormalsIntegrator, "normals")