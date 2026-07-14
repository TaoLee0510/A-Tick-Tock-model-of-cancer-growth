#pragma once

#include "config/model_config.hpp"
#include "core/types.hpp"

namespace atcg3d {

class DomainPolicy {
public:
    explicit DomainPolicy(const Model3DConfig& config)
        : bounded_(config.bounded_domain), minimum_(config.domain_min), maximum_(config.domain_max),
          thin_layer_(config.thin_layer) {}

    bool contains(Vec3i site) const noexcept {
        // Thin-layer mode constrains all cell anchors and directions to z=0.
        // z=1 remains addressable only as the upper half of a stage-0 2x2x2
        // footprint; rules never create or migrate a single-cell anchor there.
        if (thin_layer_ && site.z != 0 && site.z != 1) {
            return false;
        }
        return !bounded_ || (site.x >= minimum_.x && site.x <= maximum_.x &&
                            site.y >= minimum_.y && site.y <= maximum_.y &&
                            site.z >= minimum_.z && site.z <= maximum_.z);
    }

private:
    bool bounded_{};
    Vec3i minimum_{};
    Vec3i maximum_{};
    bool thin_layer_{};
};

}  // namespace atcg3d
