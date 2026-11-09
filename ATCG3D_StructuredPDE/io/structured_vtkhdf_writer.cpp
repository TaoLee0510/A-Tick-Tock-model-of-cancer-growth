#include "io/pde_vtkhdf_writer.hpp"
#include "model/structured_pde_model.hpp"

namespace atcg3d::continuum {
void write_pde_vtkhdf(const std::filesystem::path &path,
                      const structured_pde::StructuredPdeModel3D &model) {
    using Stage = structured_pde::StructuredStage3D;
    std::vector<GridFieldView3D> fields;
    for (auto [suffix, stage] : std::vector<std::pair<std::string, Stage>>{
             {"small", Stage::small}, {"large", Stage::large}}) {
        fields.push_back({"r_normal_" + suffix, [&model, stage](auto i) {
                              return model.r_normal(stage, i);
                          }});
        fields.push_back({"r_active_" + suffix, [&model, stage](auto i) {
                              return model.r_active(stage, i);
                          }});
        fields.push_back({"K_" + suffix, [&model, stage](auto i) {
                              return model.K(stage, i);
                          }});
        fields.push_back({"r_refractory_" + suffix, [&model, stage](auto i) {
                              return model.refractory_mass(stage, i);
                          }});
    }
    fields.push_back({"r_total", [&](auto i) {
                          return model.r_normal(Stage::small, i) +
                                 model.r_normal(Stage::large, i) +
                                 model.r_active(Stage::small, i) +
                                 model.r_active(Stage::large, i);
                      }});
    fields.push_back({"K_total", [&](auto i) {
                          return model.K(Stage::small, i) +
                                 model.K(Stage::large, i);
                      }});
    fields.push_back({"nutrient", [&](auto i) { return model.nutrient()[i]; }});
    fields.push_back({"vessel_fraction",
                      [&](auto i) { return model.vessel_fraction()[i]; }});
    fields.push_back({"occupied_fraction",
                      [&](auto i) { return model.occupied_fraction(i); }});
    if (model.angiogenesis()) {
        fields.push_back(
            {"VEGF", [&](auto i) { return model.angiogenesis()->taf()[i]; }});
        fields.push_back({"vessel_tips", [&](auto i) {
                              return model.angiogenesis()->tips()[i];
                          }});
    }
    write_grid_vtkhdf(path, model.config().continuum.grid, model.time_hours(),
                      fields);
}
} // namespace atcg3d::continuum
