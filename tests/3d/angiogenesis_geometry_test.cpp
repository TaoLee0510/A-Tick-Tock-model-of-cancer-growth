#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include "config/model_config.hpp"
#include "space/domain.hpp"
#include "vasculature/geometry.hpp"
#include "vasculature/influence_field.hpp"
#include "vasculature/surface_index.hpp"
#include "vasculature/vessel_grid.hpp"
#include "vasculature/vessel_store.hpp"

namespace {

bool close(double lhs, double rhs, double tolerance = 1e-6) {
    return std::abs(lhs - rhs) <= tolerance;
}

}  // namespace

int main() {
    using namespace atcg3d;

    VesselNodeStore3D nodes;
    VesselNodeInit3D root;
    root.position = {0, 0, 0};
    root.uid = 1;
    root.vessel_id = 11;
    root.role = VesselBranchRole::root;
    root.diameter_voxels = 2.0F;
    const VesselNodeSlot root_slot = nodes.create(root);
    assert(root_slot == 0);
    assert(nodes.position(root_slot) == (Vec3i{0, 0, 0}));
    assert(nodes.role(root_slot) == VesselBranchRole::root);
    assert(!nodes.perfused(root_slot));
    nodes.set_perfused(root_slot, true);
    assert(nodes.snapshot(root_slot).perfused);

    VesselNodeInit3D child;
    child.position = {-1, 0, 0};
    child.uid = 2;
    child.parent_uid = root.uid;
    child.parent_node_slot = root_slot;
    child.vessel_id = root.vessel_id;
    child.role = VesselBranchRole::inward;
    child.diameter_voxels = 2.0F;
    const VesselNodeSlot child_slot = nodes.create(child);
    assert(nodes.parent_node_slot(child_slot) == root_slot);
    assert(nodes.parent_uid(child_slot) == root.uid);

    VesselTipStore3D tips;
    VesselTipInit3D inward;
    inward.position = {0, 0, 0};
    inward.bias_axis = {-1, 0, 0};
    inward.target = {-10, 0, 0};
    inward.uid = 101;
    inward.vessel_id = root.vessel_id;
    inward.current_node_uid = root.uid;
    inward.current_node_slot = root_slot;
    inward.role = VesselBranchRole::inward;
    inward.status = VesselTipStatus::active;
    inward.diameter_voxels = 2.0F;
    inward.speed_voxels_per_hour = 1.5F;
    inward.max_length_voxels = 20.0F;
    inward.next_growth_time = 2.0;
    const VesselTipSlot inward_slot = tips.create(inward);
    assert(tips.current_node_slot(inward_slot) == root_slot);
    assert(tips.role(inward_slot) == VesselBranchRole::inward);
    assert(tips.status(inward_slot) == VesselTipStatus::active);
    assert(tips.consume_event_sequence(inward_slot) == 0);
    assert(tips.bump_schedule_generation(inward_slot) == 1);
    tips.set_position(inward_slot, {-1, 0, 0});
    tips.set_grown_length_voxels(inward_slot, 1.0F);
    assert(tips.snapshot(inward_slot).position == (Vec3i{-1, 0, 0}));

    bool rejected = false;
    try {
        VesselTipInit3D invalid = inward;
        invalid.uid = 102;
        invalid.speed_voxels_per_hour = 0.0F;
        tips.create(invalid);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    assert(rejected);

    const std::vector<Vec3i> thin_axis = rasterize_capsule({0, 0, 0}, {3, 0, 0}, 1.0F);
    assert(thin_axis == (std::vector<Vec3i>{{0, 0, 0}, {1, 0, 0}, {2, 0, 0}, {3, 0, 0}}));
    assert(rasterize_capsule({3, 0, 0}, {0, 0, 0}, 1.0F) == thin_axis);
    const std::vector<Vec3i> diameter_two = rasterize_capsule({0, 0, 0}, {0, 0, 0}, 2.0F);
    assert(diameter_two.size() == 7);
    assert(std::find(diameter_two.begin(), diameter_two.end(), Vec3i{1, 0, 0}) !=
           diameter_two.end());
    assert(close(segment_length({0, 0, 0}, {1, 1, 1}), std::sqrt(3.0)));

    Model3DConfig model_config;
    model_config.chunk_edge = 4;
    SparseVesselGrid3D vessel_grid(model_config.chunk_edge, DomainPolicy(model_config));
    assert(vessel_grid.add({-1, -1, -1}, VesselBranchRole::root, false, 11).placed);
    assert(vessel_grid.add({4, 0, 0}, VesselBranchRole::inward, true, 12).placed);
    assert(vessel_grid.chunk_count() == 2);
    assert(vessel_grid.occupied({-1, -1, -1}));
    assert(vessel_grid.perfused({4, 0, 0}));
    assert(vessel_grid.vessel_id({-1, -1, -1}) == 11);
    assert(vessel_grid.vessel_id({4, 0, 0}) == 12);
    const auto overlap = vessel_grid.add({4, 0, 0}, VesselBranchRole::outward, false);
    assert(overlap.placed && overlap.newly_occupied_voxels == 0);
    assert(vessel_grid.has_any({4, 0, 0}, kVesselVoxelInward));
    assert(vessel_grid.has_any({4, 0, 0}, kVesselVoxelOutward));
    assert(vessel_grid.vessel_id({4, 0, 0}) == 12);
    assert(vessel_grid.occupied_sites().size() == 2);

    Model3DConfig bounded_config;
    bounded_config.bounded_domain = true;
    bounded_config.domain_min = {-1, -1, -1};
    bounded_config.domain_max = {1, 1, 1};
    SparseVesselGrid3D bounded_grid(4, DomainPolicy(bounded_config));
    const std::array<Vec3i, 2> partly_outside{{{0, 0, 0}, {2, 0, 0}}};
    assert(!bounded_grid.add_sites(partly_outside, VesselBranchRole::root).placed);
    assert(!bounded_grid.occupied({0, 0, 0}));  // placement is preflighted atomically

    std::unordered_set<Vec3i, Vec3iHash> tumor{{0, 0, 0}};
    const auto occupied = [&tumor](Vec3i site) { return tumor.contains(site); };
    TumorSurfaceIndex3D surface;
    const std::array<Vec3i, 1> first_site{{{0, 0, 0}}};
    surface.rebuild(first_site, occupied);
    assert(surface.size() == 6);
    assert(surface.contains({{0, 0, 0}, {1, 0, 0}}));
    assert(surface.approximate_centroid() == (Vec3i{0, 0, 0}));

    tumor.insert({1, 0, 0});
    const std::array<Vec3i, 1> changed{{{1, 0, 0}}};
    surface.refresh(changed, occupied);
    assert(surface.size() == 10);
    assert(!surface.contains({{0, 0, 0}, {1, 0, 0}}));
    assert(surface.contains({{1, 0, 0}, {1, 0, 0}}));
    assert(surface.approximate_centroid() == (Vec3i{1, 0, 0}));
    const auto sample_a = surface.sample_without_replacement(4, 0.0, 77, 5);
    const auto sample_b = surface.sample_without_replacement(4, 0.0, 77, 5);
    assert(sample_a == sample_b);
    assert(sample_a.size() == 4);
    const auto complete_ranking =
        surface.sample_without_replacement(surface.size(), 0.0, 77, 5);
    assert(std::equal(sample_a.begin(), sample_a.end(), complete_ranking.begin()));
    assert(surface.sample_without_replacement(4, 100.0, 77, 5).size() == 1);

    // A stable-hash top-1 sample should not privilege a face because of the
    // unordered_set iteration order. Use a deliberately generous tolerance
    // to catch fixed-axis bias without making this a flaky RNG quality test.
    const auto ordered_faces = surface.faces();
    std::vector<std::size_t> selection_counts(ordered_faces.size(), 0U);
    constexpr std::size_t kSamplingTrials = 10000U;
    for (std::size_t trial = 0; trial < kSamplingTrials; ++trial) {
        const auto selected = surface.sample_without_replacement(1, 0.0, 77, trial);
        assert(selected.size() == 1);
        const auto position =
            std::find(ordered_faces.begin(), ordered_faces.end(), selected.front());
        assert(position != ordered_faces.end());
        ++selection_counts[static_cast<std::size_t>(position - ordered_faces.begin())];
    }
    const double expected =
        static_cast<double>(kSamplingTrials) / static_cast<double>(ordered_faces.size());
    for (const std::size_t observed : selection_counts) {
        assert(static_cast<double>(observed) > expected * 0.75);
        assert(static_cast<double>(observed) < expected * 1.25);
    }

    tumor.erase({1, 0, 0});
    surface.refresh(changed, occupied);
    assert(surface.size() == 6);
    assert(surface.contains({{0, 0, 0}, {1, 0, 0}}));
    assert(surface.approximate_centroid() == (Vec3i{0, 0, 0}));

    TumorSurfaceIndex3D rebuilt_surface;
    rebuilt_surface.rebuild(first_site, occupied);
    assert(rebuilt_surface.faces() == surface.faces());
    assert(rebuilt_surface.approximate_centroid() == surface.approximate_centroid());

    // Lesion-specific sampling must build axial extrema from the selected
    // component only. A remote component is neither a candidate nor an
    // occluder for the requested lesion.
    std::unordered_set<Vec3i, Vec3iHash> two_components{{0, 0, 0}, {10, 0, 0}};
    const std::vector<Vec3i> two_component_sites{{0, 0, 0}, {10, 0, 0}};
    TumorSurfaceIndex3D component_surface;
    component_surface.rebuild(two_component_sites,
        [&two_components](Vec3i site) { return two_components.contains(site); });
    const auto first_component_faces =
        component_surface.sample_external_subset_without_replacement(
            component_surface.size(), 0.0, 7, 0,
            [](const ExposedFace3D& face) { return face.inside.x < 5; });
    assert(first_component_faces.size() == 6U);
    assert(std::all_of(first_component_faces.begin(), first_component_faces.end(),
                       [](const ExposedFace3D& face) {
                           return face.inside == Vec3i{0, 0, 0};
                       }));

    // A one-voxel-thick 5x5x5 shell exposes both its outside and its enclosed
    // 3x3x3 cavity. Axis-visible sampling must retain exactly the 6*5*5 outer
    // faces and must never select an inward-facing cavity wall.
    std::unordered_set<Vec3i, Vec3iHash> hollow_tumor;
    std::vector<Vec3i> hollow_sites;
    for (std::int32_t z = -2; z <= 2; ++z) {
        for (std::int32_t y = -2; y <= 2; ++y) {
            for (std::int32_t x = -2; x <= 2; ++x) {
                if (std::abs(x) != 2 && std::abs(y) != 2 && std::abs(z) != 2) continue;
                const Vec3i site{x, y, z};
                hollow_tumor.insert(site);
                hollow_sites.push_back(site);
            }
        }
    }
    const auto hollow_occupied = [&hollow_tumor](Vec3i site) {
        return hollow_tumor.contains(site);
    };
    TumorSurfaceIndex3D hollow_surface;
    hollow_surface.rebuild(hollow_sites, hollow_occupied);
    assert(hollow_surface.size() == 204U);  // 150 outer + 54 cavity faces
    const auto is_outer_shell_face = [](const ExposedFace3D& face) {
        const Vec3i normal = face.outward_normal;
        if (normal.x != 0) return face.inside.x == 2 * normal.x;
        if (normal.y != 0) return face.inside.y == 2 * normal.y;
        return face.inside.z == 2 * normal.z;
    };
    const auto all_hollow_external = hollow_surface.sample_external_without_replacement(
        hollow_surface.size(), 0.0, 90210, 4);
    assert(all_hollow_external.size() == 150U);
    assert(std::all_of(all_hollow_external.begin(), all_hollow_external.end(),
                       is_outer_shell_face));
    for (std::uint64_t sequence = 0; sequence < 128U; ++sequence) {
        const auto sampled =
            hollow_surface.sample_external_without_replacement(16, 0.0, 90210, sequence);
        assert(sampled.size() == 16U);
        assert(std::all_of(sampled.begin(), sampled.end(), is_outer_shell_face));
    }

    // Every face of a solid cube is external. The external sampler must match
    // the original stable-hash sample exactly and remain symmetric across the
    // six axis directions over independent event sequences.
    std::unordered_set<Vec3i, Vec3iHash> solid_tumor;
    std::vector<Vec3i> solid_sites;
    for (std::int32_t z = -1; z <= 1; ++z) {
        for (std::int32_t y = -1; y <= 1; ++y) {
            for (std::int32_t x = -1; x <= 1; ++x) {
                const Vec3i site{x, y, z};
                solid_tumor.insert(site);
                solid_sites.push_back(site);
            }
        }
    }
    const auto solid_occupied = [&solid_tumor](Vec3i site) {
        return solid_tumor.contains(site);
    };
    TumorSurfaceIndex3D solid_surface;
    solid_surface.rebuild(solid_sites, solid_occupied);
    assert(solid_surface.size() == 54U);
    const auto solid_external_a =
        solid_surface.sample_external_without_replacement(8, 0.0, 314159, 12);
    const auto solid_external_b =
        solid_surface.sample_external_without_replacement(8, 0.0, 314159, 12);
    assert(solid_external_a == solid_external_b);
    assert(solid_external_a ==
           solid_surface.sample_without_replacement(8, 0.0, 314159, 12));

    const auto normal_bucket = [](Vec3i normal) -> std::size_t {
        if (normal.x < 0) return 0U;
        if (normal.x > 0) return 1U;
        if (normal.y < 0) return 2U;
        if (normal.y > 0) return 3U;
        if (normal.z < 0) return 4U;
        return 5U;
    };
    std::array<std::size_t, 6> direction_counts{};
    constexpr std::size_t kExternalSamplingTrials = 12000U;
    for (std::uint64_t sequence = 0; sequence < kExternalSamplingTrials; ++sequence) {
        const auto sampled =
            solid_surface.sample_external_without_replacement(1, 0.0, 314159, sequence);
        assert(sampled.size() == 1U);
        ++direction_counts[normal_bucket(sampled.front().outward_normal)];
    }
    const double expected_direction_count =
        static_cast<double>(kExternalSamplingTrials) /
        static_cast<double>(direction_counts.size());
    for (const std::size_t observed : direction_counts) {
        assert(static_cast<double>(observed) > expected_direction_count * 0.85);
        assert(static_cast<double>(observed) < expected_direction_count * 1.15);
    }

    VascularInfluenceField3D influence(4, 2.0F, 0.8F);
    influence.add_source({-1, 0, 0});
    assert(close(influence.relief({-1, 0, 0}), 0.8));
    assert(close(influence.relief({0, 0, 0}), 0.4));
    assert(close(influence.relief({1, 0, 0}), 0.0));
    assert(close(influence.effective_density(20.0, {-1, 0, 0}), 4.0));
    assert(close(influence.local_capacity_multiplier({-1, 0, 0}), 5.0));
    const float before_overlap = influence.relief({0, 0, 0});
    influence.add_source({-1, 1, 0});
    assert(influence.relief({0, 0, 0}) >= before_overlap);
    assert(influence.relief({0, 0, 0}) <= 0.8F);
    assert(influence.block_count() >= 2);  // the source neighborhood crosses x=-1/0 blocks

    return 0;
}
