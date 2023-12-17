#include <lightwave.hpp>

#include "fresnel.hpp"
#include "microfacet.hpp"

namespace lightwave {

    struct DiffuseLobe {
        Color color;

        BsdfEval evaluate(const Vector& wo, const Vector& wi) const {
            // hints:
            // * copy your diffuse bsdf evaluate here
            // * you do not need to query a texture, the albedo is given by `color`
            return { color * max(0.0f, Frame::cosTheta(wi)) / Pi };
        }

        BsdfSample sample(const Vector& wo, Sampler& rng) const {
            // hints:
            // * copy your diffuse bsdf evaluate here
            // * you do not need to query a texture, the albedo is given by `color`
            Point2 sampledPoint = rng.next2D(); 
            Vector sampledVector = squareToCosineHemisphere(sampledPoint);
            return BsdfSample(sampledVector.normalized(), color);
        }
    };

    struct MetallicLobe {
        float alpha;
        Color color;

        BsdfEval evaluate(const Vector& wo, const Vector& wi) const {
            // hints:
            // * copy your roughconductor bsdf evaluate here
            // * you do not need to query textures
            //   * the reflectance is given by `color'
            //   * the variable `alpha' is already provided for you
            Vector wm = (wi + wo).normalized();
            float D = microfacet::evaluateGGX(alpha, wm);
            float G1wi = microfacet::smithG1(alpha, wm, wi);
            float G1wo = microfacet::smithG1(alpha, wm, wo);
            float denominator = 4 * Frame::absCosTheta(wi) * Frame::absCosTheta(wo);
            return BsdfEval(color * D * G1wi * G1wo / denominator * Frame::cosTheta(wi));
        }

        BsdfSample sample(const Vector& wo, Sampler& rng) const {
            // hints:
            // * copy your roughconductor bsdf sample here
            // * you do not need to query textures
            //   * the reflectance is given by `color'
            //   * the variable `alpha' is already provided for you
            Vector m = microfacet::sampleGGXVNDF(alpha, wo, rng.next2D()); // Sample microfacet normal m
            Vector wi = reflect(wo, m); // Reflect wo about m to get outgoing direction wi
            Vector wm = (wi + wo).normalized();
            float G1wi = microfacet::smithG1(alpha, wm, wi);
            return BsdfSample(wi, color * G1wi);
        }
    };

    class Principled : public Bsdf {
        ref<Texture> m_baseColor;
        ref<Texture> m_roughness;
        ref<Texture> m_metallic;
        ref<Texture> m_specular;

        struct Combination {
            float diffuseSelectionProb;
            DiffuseLobe diffuse;
            MetallicLobe metallic;
        };

        Combination combine(const Point2& uv, const Vector& wo) const {
            const auto baseColor = m_baseColor->evaluate(uv);
            const auto alpha = std::max(float(1e-3), sqr(m_roughness->scalar(uv)));
            const auto specular = m_specular->scalar(uv);
            const auto metallic = m_metallic->scalar(uv);
            const auto F =
                specular * schlick((1 - metallic) * 0.08f, Frame::cosTheta(wo));

            const DiffuseLobe diffuseLobe = {
                .color = (1 - F) * (1 - metallic) * baseColor,
            };
            const MetallicLobe metallicLobe = {
                .alpha = alpha,
                .color = F * Color(1) + (1 - F) * metallic * baseColor,
            };

            const auto diffuseAlbedo = diffuseLobe.color.mean();
            const auto totalAlbedo =
                diffuseLobe.color.mean() + metallicLobe.color.mean();
            return {
                .diffuseSelectionProb =
                    totalAlbedo > 0 ? diffuseAlbedo / totalAlbedo : 1.0f,
                .diffuse = diffuseLobe,
                .metallic = metallicLobe,
            };
        }

        public:
        Principled(const Properties& properties) {
            m_baseColor = properties.get<Texture>("baseColor");
            m_roughness = properties.get<Texture>("roughness");
            m_metallic = properties.get<Texture>("metallic");
            m_specular = properties.get<Texture>("specular");
        }

        BsdfEval evaluate(const Point2& uv, const Vector& wo,
                          const Vector& wi) const override {
            const auto combination = combine(uv, wo);

            // hint: evaluate `combination.diffuse` and `combination.metallic` and
            // combine their results
            BsdfEval diffuseEval = combination.diffuse.evaluate(wo, wi);
            BsdfEval metallicEval = combination.metallic.evaluate(wo, wi);
            return BsdfEval(diffuseEval.value + metallicEval.value);
        }

        BsdfSample sample(const Point2& uv, const Vector& wo,
                          Sampler& rng) const override {
            const auto combination = combine(uv, wo);

            // hint: sample either `combination.diffuse` (probability
            // `combination.diffuseSelectionProb`) or `combination.metallic`
            if (rng.next() < combination.diffuseSelectionProb) {
                // Sample diffuse lobe
                BsdfSample sample = combination.diffuse.sample(wo, rng);
                return BsdfSample(sample.wi, sample.weight / combination.diffuseSelectionProb);
            } else {
                // Sample metallic lobe
                BsdfSample sample = combination.metallic.sample(wo, rng);
                return BsdfSample(sample.wi, sample.weight / (1.0f - combination.diffuseSelectionProb));
            }
        }

        std::string toString() const override {
            return tfm::format("Principled[\n"
                               "  baseColor = %s,\n"
                               "  roughness = %s,\n"
                               "  metallic  = %s,\n"
                               "  specular  = %s,\n"
                               "]",
                               indent(m_baseColor), indent(m_roughness),
                               indent(m_metallic), indent(m_specular));
        }
    };
} // namespace lightwave

REGISTER_BSDF(Principled, "principled")