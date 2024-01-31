#include <lightwave/camera.hpp>
#include <lightwave/sampler.hpp>

namespace lightwave
{
    CameraSample Camera::sample(const Point2i &pixel, Sampler &rng) const
    {
        // begin by sampling a random position within the pixel
        const auto pixelPlusRandomOffset = Vector2(pixel.cast<float>()) + Vector2(rng.next2D());
        // normalize by image resolution to end up with value in range [-1,-1] to [+1,+1]
        const auto normalized = 2 * pixelPlusRandomOffset / m_resolution.cast<float>() - Vector2(1);
        // generate the sample using the normalized sample function
        const auto s = sample(normalized, rng);
        assert_normalized(s.ray.direction, {
            logger(EError, "your Camera::sample implementation returned a non-normalized direction");
        });
        return s;
    }

    void buildBokehShapes(Bokeh &bokehConfig)
    {
        int nVert = bokehConfig.blades;
        if (nVert < 3) // No polygons with fewer than 3 vertices
            return;

        bokehConfig.step = 2 * Pi / nVert; // Calculate angular step based on the current number of vertices
        bokehConfig.corners = std::unique_ptr<Point2[]>(new Point2[nVert]);
        float cur = Pi2; // Start at (0, 1) with all polygons

        // Populate with corners of the current polygon
        for (int i = 0; i < nVert; i++)
        {
            bokehConfig.corners[i] = Point2(cos(cur), sin(cur));
            cur += bokehConfig.step;
        }
    }
} // namespace lightwave
