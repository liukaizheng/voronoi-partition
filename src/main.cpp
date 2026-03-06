#include <iostream>
#include <CLI/CLI.hpp>
#include <fstream>
#include <boost/functional/hash.hpp>
#include <ranges>
#include <deque>
#include <algorithm>
#include <limits>
#include <cmath>
#include <cstdint>
#include <iomanip>

#include <Eigen/Dense>
#include <gpf/mesh.hpp>
#include <unordered_set>
#include <vector>
#include <voro++.hh>

namespace ranges = std::ranges;
namespace views = ranges::views;

auto read_colored_mesh(const std::string& file_name) {
    std::ifstream input_file(file_name, std::ios::in);
    int n_v, n_f, n_e;
    std::string line;
    std::vector<std::array<double, 3>> points;
    std::vector<std::vector<std::size_t>> faces;
    getline(input_file, line); // file format
    {
        line.clear();
        getline(input_file, line);
        std::stringstream line_stream;
        line_stream.str(line);
        line_stream >> n_v >> n_f >> n_e;
    }

    std::vector<std::array<std::size_t, 3>> face_colors(n_f);
    std::array<double, 3> pos;
    for(int i = 0; i < n_v; ++i) {
        line.clear();
        getline(input_file, line);
        std::stringstream line_stream;
        line_stream.str(line);
        line_stream >> pos[0] >> pos[1] >> pos[2];
        points.emplace_back(pos);
    }

    std::size_t indices[3];
    std::size_t color[3];
    int f;
    for(int i = 0; i < n_f; ++i) {
        line.clear();
        getline(input_file, line);
        std::stringstream line_stream;
        line_stream.str(line);
        line_stream >> f;
        for(int j = 0; j < f; ++j) {
        line_stream >> indices[j];
        }
        line_stream >> color[0] >> color[1] >> color[2];
        face_colors[i] = {color[0], color[1], color[2]};
        faces.emplace_back(std::vector{indices[0], indices[1], indices[2]});
    }
    input_file.close();
    auto hash = [](const std::array<std::size_t, 3>& key) {
        return boost::hash_value(key);
    };
    std::unordered_map<std::array<std::size_t, 3>, std::vector<std::size_t>, decltype(hash)> color_to_face_map(faces.size(), hash);
    for (std::size_t i = 0; i < faces.size(); ++i) {
        color_to_face_map[face_colors[i]].push_back(i);
    }
    return std::make_tuple(
        std::move(points),
        std::move(faces),
        std::move(color_to_face_map)
            | std::views::elements<1>
            | views::transform([](auto&& faces) {
                return faces | views::transform([](const auto fid) {
                    return gpf::FaceId{fid}; }) | ranges::to<std::vector>();
                })
            | ranges::to<std::vector>()
    );
}

auto extract_color_boundaries(
    const std::vector<std::array<double, 3>>& points,
    const std::vector<std::vector<std::size_t>>& faces,
    std::vector<std::vector<gpf::FaceId>> color_face_groups
)  {

    struct EdgeProp {
        bool is_boundary = false;
    };
    struct FaceProp {
        std::size_t color_index;
        std::size_t region_index = gpf::kInvalidIndex;
    };

    using Mesh = gpf::ManifoldMesh<gpf::Empty, gpf::Empty, EdgeProp, FaceProp>;
    auto mesh = Mesh::new_in(faces);
    for (std::size_t i = 0; i < color_face_groups.size(); i++) {
        for (const auto fid : color_face_groups[i]) {
            mesh.face_prop(fid).color_index = i;
        }
    }
    for (auto edge : mesh.edges()) {
        auto ha = edge.halfedge();
        auto fa = ha.face();
        auto fb = ha.twin().face();
        if (fa.prop().color_index != fb.prop().color_index) {
            edge.prop().is_boundary = true;
        }
    }

    std::vector<std::size_t> region_colors;
    for (auto face : mesh.faces()) {
        if (face.prop().region_index != gpf::kInvalidIndex) {
            continue;
        }
        const auto region_index = region_colors.size();
        face.prop().region_index = region_index;
        region_colors.push_back(face.prop().color_index);
        std::deque<gpf::FaceId> queue{face.id};
        while (!queue.empty()) {
            auto curr_fid = queue.front();
            queue.pop_front();
            for (auto he: mesh.face(curr_fid).halfedges()) {
                if (he.edge().prop().is_boundary) {
                    continue;
                }
                auto adj_face = he.twin().face();
                auto& adj_face_prop = adj_face.prop();
                if (adj_face_prop.region_index != gpf::kInvalidIndex) {
                    continue;
                }
                adj_face_prop.region_index = region_index;
                queue.push_back(adj_face.id);
            }
        }
    }
    return std::make_pair(std::move(mesh), std::move(region_colors));
}


template <class Mesh>
auto compute_edge_offset_points(
    const std::vector<std::array<double, 3>>& points,
    const Mesh& mesh,
    const std::size_t n_regions,
    double offset = 0.0001
) {
    using Vector3 = Eigen::Vector3d;
    std::vector<std::array<double, 3>> seed_points;
    std::vector<std::size_t> point_group_indices;

    auto compute_face_normal = [&points](auto face) -> Vector3 {
        std::array<const double*, 3> pts;
        int i = 0;
        for (auto he : face.halfedges()) {
            pts[i++] = points[he.to().id.idx].data();
        }
        auto p1 = Vector3::Map(pts[0]);
        auto p2 = Vector3::Map(pts[1]);
        auto p3 = Vector3::Map(pts[2]);
        return (p2 - p1).cross(p3 - p1).normalized();
    };

    auto push_point = [&seed_points, &point_group_indices](const Vector3& pt, std::size_t region_idx) {
        seed_points.emplace_back();
        point_group_indices.emplace_back(region_idx);
        Vector3::Map(seed_points.back().data()) = pt;
    };

    std::vector<bool> visited(mesh.n_vertices_capacity(), false);
    constexpr double step = 0.5;
    for (auto edge : mesh.edges()) {
        auto ha = edge.halfedge();
        auto hb = ha.twin();
        auto fa = ha.face();
        auto fb = hb.face();

        auto region_a = fa.prop().region_index;
        auto region_b = fb.prop().region_index;


        if (region_a != region_b) {
            auto va = hb.to().id.idx;
            auto vb = ha.to().id.idx;
            visited[va] = true;
            visited[vb] = true;

            auto pa = Vector3::Map(points[va].data());
            auto pb = Vector3::Map(points[vb].data());
            auto mid_pt = ((pa + pb) * 0.5).eval();
            auto edge_dir = (pb - pa).eval();

            auto normal_a = compute_face_normal(fa);
            auto dir_a = normal_a.cross(edge_dir).normalized().eval();
            // push_point(mid_pt + dir_a * offset, region_a);

            auto normal_b = compute_face_normal(fb);
            auto dir_b = normal_b.cross(-edge_dir).normalized().eval();
            // push_point(mid_pt + dir_b * offset, region_b);
            const std::size_t n_samples = std::max(2ul, static_cast<std::size_t>(std::round(edge_dir.norm() / step)));
            {
                auto pt = (pa + edge_dir * 0.01).eval();
                push_point(pt + dir_a * offset, region_a);
                push_point(pt + dir_b * offset, region_b);
                pt = (pa + edge_dir * 0.99).eval();
                push_point(pt + dir_a * offset, region_a);
                push_point(pt + dir_b * offset, region_b);
            }
            for (int i = 0; i < n_samples; ++i) {
                double t = (i + 0.5) / n_samples;
                auto pt = (pa + t * edge_dir).eval();
                push_point(pt + dir_a * offset, region_a);
                push_point(pt + dir_b * offset, region_b);
            }
        }

    }
    for (auto v : mesh.vertices()) {
        if (visited[v.id.idx])  {
            continue;
        }
        visited[v.id.idx] = true;

        auto face = v.halfedge().face();
        seed_points.push_back(points[v.id.idx]);
        point_group_indices.push_back(face.prop().region_index);
    }

    return std::make_pair(std::move(seed_points), std::move(point_group_indices));
}

struct VoronoiCell {
    int id;
    double x, y, z;
    double volume;
    std::vector<double> vertices;
    std::vector<int> face_vertices;
    std::vector<int> neighbors;
};

auto compute_voronoi(
    const std::vector<std::array<double, 3>>& seed_points
) {
    // Compute bounding box with padding
    double min_x = std::numeric_limits<double>::max();
    double min_y = min_x, min_z = min_x;
    double max_x = std::numeric_limits<double>::lowest();
    double max_y = max_x, max_z = max_x;

    for (const auto& p : seed_points) {
        min_x = std::min(min_x, p[0]); max_x = std::max(max_x, p[0]);
        min_y = std::min(min_y, p[1]); max_y = std::max(max_y, p[1]);
        min_z = std::min(min_z, p[2]); max_z = std::max(max_z, p[2]);
    }

    // Add padding so boundary cells are well-formed
    // double pad_x = (max_x - min_x) * 0.0 + 1e-6;
    // double pad_y = (max_y - min_y) * 0.0 + 1e-6;
    // double pad_z = (max_z - min_z) * 0.0 + 1e-6;
    // min_x -= pad_x; max_x += pad_x;
    // min_y -= pad_y; max_y += pad_y;
    // min_z -= pad_z; max_z += pad_z;

    // Choose grid resolution based on particle count
    int n = static_cast<int>(seed_points.size());
    int n_blocks = std::max(1, static_cast<int>(std::cbrt(n / 5.0)));

    voro::container con(
        min_x, max_x, min_y, max_y, min_z, max_z,
        n_blocks, n_blocks, n_blocks,
        false, false, false, 8
    );

    for (int i = 0; i < n; ++i) {
        con.put(i, seed_points[i][0], seed_points[i][1], seed_points[i][2]);
    }

    std::vector<VoronoiCell> cells;
    cells.reserve(n);

    voro::voronoicell_neighbor c;
    voro::c_loop_all loop(con);
    if (loop.start()) do {
        if (con.compute_cell(c, loop)) {
            VoronoiCell cell;
            cell.id = loop.pid();
            loop.pos(cell.x, cell.y, cell.z);
            cell.volume = c.volume();
            c.vertices(cell.x, cell.y, cell.z, cell.vertices);
            c.face_vertices(cell.face_vertices);
            c.neighbors(cell.neighbors);
            cells.push_back(std::move(cell));
        }
    } while (loop.inc());

    // Sort by id so the output order matches the input seed order
    std::sort(cells.begin(), cells.end(),
        [](const VoronoiCell& a, const VoronoiCell& b) { return a.id < b.id; });

    return cells;
}

struct InterfaceMesh {
    std::vector<std::array<double, 3>> vertices;
    std::vector<std::vector<std::size_t>> faces;
    std::vector<std::array<std::size_t, 2>> face_groups;
};

auto extract_interface_mesh(
    const std::vector<VoronoiCell>& cells,
    const std::vector<std::size_t>& point_group_indices,
    double tolerance = 1e-3
) -> InterfaceMesh {
    InterfaceMesh result;

    double inv_tol = 1.0 / tolerance;
    double tol_sq = tolerance * tolerance;

    auto hasher = [](const std::array<int64_t, 3>& key) -> std::size_t {
        return boost::hash_value(key);
    };

    std::unordered_map<std::array<int64_t, 3>, std::size_t, decltype(hasher)> vertex_map(0, hasher);

    auto get_or_insert = [&](double x, double y, double z) -> std::size_t {
        int64_t ix = static_cast<int64_t>(std::floor(x * inv_tol));
        int64_t iy = static_cast<int64_t>(std::floor(y * inv_tol));
        int64_t iz = static_cast<int64_t>(std::floor(z * inv_tol));

        // Check current cell and all 26 neighbors to handle boundary cases
        for (int dx = -1; dx <= 1; ++dx) {
            for (int dy = -1; dy <= 1; ++dy) {
                for (int dz = -1; dz <= 1; ++dz) {
                    auto it = vertex_map.find({ix + dx, iy + dy, iz + dz});
                    if (it != vertex_map.end()) {
                        const auto& existing = result.vertices[it->second];
                        double ex = x - existing[0];
                        double ey = y - existing[1];
                        double ez = z - existing[2];
                        if (ex*ex + ey*ey + ez*ez < tol_sq) {
                            return it->second;
                        }
                    }
                }
            }
        }

        std::size_t idx = result.vertices.size();
        vertex_map.emplace(std::array<int64_t, 3>{ix, iy, iz}, idx);
        result.vertices.push_back({x, y, z});
        return idx;
    };

    for (const auto& cell : cells) {
        std::size_t my_group = point_group_indices[cell.id];

        for (std::size_t i = 0, j = 0; i < cell.neighbors.size(); ++i) {
            int n_verts = cell.face_vertices[j];
            int neighbor_id = cell.neighbors[i];

            bool is_boundary = neighbor_id < 0;
            bool is_interface = neighbor_id >= 0
                && cell.id < neighbor_id
                && static_cast<std::size_t>(neighbor_id) < point_group_indices.size()
                && my_group != point_group_indices[neighbor_id];

            if (is_boundary || is_interface) {
                std::vector<std::size_t> face;
                face.reserve(n_verts);
                for (int k = 0; k < n_verts; ++k) {
                    int vi = cell.face_vertices[j + 1 + k];
                    auto idx = get_or_insert(
                        cell.vertices[3 * vi],
                        cell.vertices[3 * vi + 1],
                        cell.vertices[3 * vi + 2]
                    );
                    if (face.empty() || face.back() != idx) {
                        face.push_back(idx);
                    }
                }
                // Remove wrap-around duplicate
                if (face.size() > 1 && face.front() == face.back()) {
                    face.pop_back();
                }
                auto face_set = face | ranges::to<std::unordered_set>();
                if (face_set.size() != face.size()) {
                    const auto a = 2;
                }
                if (face.size() >= 3) {
                    result.faces.push_back(std::move(face));
                    if (is_boundary) {
                        result.face_groups.push_back(
                            {gpf::kInvalidIndex, my_group});
                    } else {
                        result.face_groups.push_back(
                            {point_group_indices[neighbor_id], my_group});
                    }
                }
            }

            j += n_verts + 1;
        }
    }

    return result;
}

void write_seed_points(const std::string& path,
                       const std::vector<std::array<double, 3>>& seed_points,
                       const std::vector<std::size_t>& point_group_indices) {
    std::ofstream out(path);
    out << "OFF\n";
    out << seed_points.size() << " 0 0\n";
    out << std::setprecision(17);
    for (const auto& p : seed_points) {
        out << p[0] << " " << p[1] << " " << p[2] << "\n";
    }
}

void write_off(const std::string& path, const InterfaceMesh& mesh) {
    std::ofstream out(path);
    out << "OFF\n";
    out << mesh.vertices.size() << " " << mesh.faces.size() << " 0\n";
    out << std::setprecision(17);
    for (const auto& v : mesh.vertices) {
        out << v[0] << " " << v[1] << " " << v[2] << "\n";
    }
    for (const auto& f : mesh.faces) {
        out << f.size();
        for (auto idx : f) out << " " << idx;
        out << "\n";
    }
}

int main(int argc, char** argv) {
    CLI::App app { "multi-label-partition" };
    std::string mesh_path;
    std::string output_path = "interface.off";
    app.add_option("-m,--mesh", mesh_path, "Path to mesh file")->required();
    app.add_option("-o,--output", output_path, "Path to output interface mesh");
    std::string seed_path = "seed_points.off";
    app.add_option("-s,--seeds", seed_path, "Path to output seed points");
    CLI11_PARSE(app, argc, argv);

    std::vector<std::array<std::size_t, 3>> colors;
    auto [points, faces, label_face_groups] = read_colored_mesh(mesh_path);

    auto [mesh, region_colors] = extract_color_boundaries(points, faces, label_face_groups);
    auto [seed_points,  point_group_indices] = compute_edge_offset_points(points, mesh, region_colors.size());
    write_seed_points(seed_path, seed_points, point_group_indices);

    auto cells = compute_voronoi(seed_points);

    auto interface_mesh = extract_interface_mesh(cells, point_group_indices);
    write_off(output_path, interface_mesh);

    std::cout << "Computed " << cells.size() << " Voronoi cells from "
              << seed_points.size() << " seed points\n";
    std::cout << "Interface mesh: " << interface_mesh.vertices.size() << " vertices, "
              << interface_mesh.faces.size() << " faces\n";
    std::cout << "Written to: " << output_path << "\n";
}
