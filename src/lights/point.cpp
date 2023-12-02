#include <lightwave.hpp>

namespace lightwave {

class PointLight final : public Light {
public:
    PointLight(const Properties &properties) {
    }

    DirectLightSample sampleDirect(const Point &origin,
                                   Sampler &rng) const override {
        NOT_IMPLEMENTED
    }

    bool canBeIntersected() const override { return false; }

    std::string toString() const override {
        return tfm::format("PointLight[\n"
                           "]");
    }
};

} // namespace lightwave

REGISTER_LIGHT(PointLight, "point")
