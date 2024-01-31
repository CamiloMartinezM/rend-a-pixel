#include <OpenImageDenoise/oidn.hpp>
#include <lightwave.hpp>

namespace lightwave
{
    class Denoising : public Postprocess
    {
      protected:
        /// @brief Holds the optional normals AOV.
        ref<Image> m_normals;
        /// @brief Holds the optional albedo AOV.
        ref<Image> m_albedo;

      public:
        Denoising(const Properties &properties) : Postprocess(properties)
        {
            m_normals = properties.get<Image>("normals", nullptr);
            m_albedo = properties.get<Image>("albedo", nullptr);
        }

        void execute() override
        {
            const Point2i resolution = m_input->resolution();
            const int width = resolution.x();
            const int height = resolution.y();
            m_output->initialize(resolution);
            if (!m_input || !m_output)
                throw std::runtime_error("Input or output image not set for Denoising");

            // Initialize OIDN device
            oidn::DeviceRef device = oidn::newDevice(oidn::DeviceType::CPU);
            device.commit();

            // Create and set up the denoising filter
            oidn::FilterRef filter = device.newFilter("RT");
            filter.setImage("color", m_input->data(), oidn::Format::Float3, width, height);

            // Add normals and albedo to the filter as AOVs, if they are available
            // If only normals is provided, it throws the Error: unsupported combination of input features
            // https://github.com/OpenImageDenoise/oidn/issues/9
            if (m_normals && m_albedo)
            {
                filter.setImage("normal", m_normals->data(), oidn::Format::Float3, width, height);
                filter.setImage("albedo", m_albedo->data(), oidn::Format::Float3, width, height);
            }

            // Set up the output
            filter.setImage("output", m_output->data(), oidn::Format::Float3, width, height);
            filter.set("hdr", true);
            filter.setProgressMonitorFunction([](void *userPtr, double n) -> bool {
                logger(EDebug, "Denoising image...\n");
                return true;
            });

            // Execute the denoiser
            filter.commit();
            filter.execute();

            // Check for errors if any
            const char *errorMessage;
            if (device.getError(errorMessage) != oidn::Error::None)
                throw std::runtime_error(errorMessage);

            logger(EDebug, "Saving denoised image...\n");
            m_output->save(); // Save the output image

            // Show the image on Tev
            Streaming stream { *m_output };
            stream.update();
        }

        /// @brief Returns a textual representation of the Denoising AOVs.
        std::string toString() const override
        {
            return tfm::format("Denoising[\n"
                               "  m_input = %s,\n"
                               "  m_output = %s,\n"
                               "  m_normals = %s,\n"
                               "  m_albedo = %s,\n"
                               "]",
                               indent(m_input), indent(m_output), indent(m_normals), indent(m_albedo));
        }
    };
} // namespace lightwave

// Register the Denoising class for the "denoising" postprocess type.
REGISTER_POSTPROCESS(Denoising, "denoising")